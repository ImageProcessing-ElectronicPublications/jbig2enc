// Copyright 2006 Google Inc. All Rights Reserved.
// Author: agl@imperialviolet.org (Adam Langley)
//
// Copyright (C) 2006 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <vector>

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include <leptonica/allheaders.h>
#if (LIBLEPT_MAJOR_VERSION == 1 && LIBLEPT_MINOR_VERSION >= 83) || LIBLEPT_MAJOR_VERSION > 1
#include "leptonica/pix_internal.h"
#endif

#include "jbig2enc.h"

#if defined(WIN32)
#define WINBINARY O_BINARY
#else
#define WINBINARY 0
#endif

#define JBIG2_THRESHOLD_MIN 0.4f
#define JBIG2_THRESHOLD_MAX 0.97f
#define JBIG2_THRESHOLD_DEF 0.92f
#define JBIG2_WEIGHT_MIN 0.1f
#define JBIG2_WEIGHT_MAX 0.9f
#define JBIG2_WEIGHT_DEF 0.5f
#define BW_THRESHOLD_MIN 0
#define BW_THRESHOLD_MAX 255
#define BW_LOCAL_THRESHOLD_DEF 200
#define BW_GLOBAL_THRESHOLD_DEF 128

static void
usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [options] <input filenames...>\n", argv0);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -d --duplicate-line-removal: use TPGD in generic region coder\n");
  fprintf(stderr, "  -p --pdf: produce PDF ready data\n");
  fprintf(stderr, "  -s --symbol-mode: use text region, not generic coder\n");
  fprintf(stderr, "  -t <threshold>: set classification threshold for symbol coder (def: %0.2f)\n", JBIG2_THRESHOLD_DEF);
  fprintf(stderr, "  -w <weight>: set classification weight for symbol coder (def: %0.2f)\n", JBIG2_WEIGHT_DEF);
  fprintf(stderr, "  -P <number> --pages-per-dict <number>: pages per dictionary (default 15)\n");
  fprintf(stderr, "  -T <bw threshold>: set 1 bpp threshold (def: %d)\n", BW_LOCAL_THRESHOLD_DEF);
  fprintf(stderr, "  -G --global: use global BW threshold on 8 bpp images;\n"
                  "               the default is to use local (adaptive) thresholding\n");
  fprintf(stderr, "  -r --refine: use refinement (requires -s: lossless)\n");
  fprintf(stderr, "  -O <outfile>: dump thresholded image as PNG\n");
  fprintf(stderr, "  -2: upsample 2x before thresholding\n");
  fprintf(stderr, "  -4: upsample 4x before thresholding\n");
  fprintf(stderr, "  -S: remove images from mixed input and save separately\n");
  fprintf(stderr, "  -j --jpeg-output: write images from mixed input as JPEG\n");
  fprintf(stderr, "  -a --auto-thresh: use automatic thresholding in symbol encoder\n");
  fprintf(stderr, "  -D --dpi: force dpi\n");
  fprintf(stderr, "  --no-hash: disables use of hash function for automatic thresholding\n");
  fprintf(stderr, "  -V --version: version info\n");
  fprintf(stderr, "  -v: be verbose\n");
}

static bool verbose = false;


static void
pixInfo(PIX *pix, const char *msg)
{
  if (msg != NULL) fprintf(stderr, "%s ", msg);
  if (pix == NULL)
  {
    fprintf(stderr, "NULL pointer!\n");
    return;
  }
  fprintf(stderr, "%u x %u (%d bits) %udpi x %udpi, refcount = %u\n",
          pix->w, pix->h, pix->d, pix->xres, pix->yres, pix->refcount);
}

#ifdef WIN32
// -----------------------------------------------------------------------------
// Windows, sadly, lacks asprintf
// -----------------------------------------------------------------------------
#include <stdarg.h>
int
asprintf(char **strp, const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);

  const int required = vsnprintf(NULL, 0, fmt, va);
  char *const buffer = (char *) malloc(required + 1);
  const int ret = vsnprintf(buffer, required + 1, fmt, va);
  *strp = buffer;

  va_end(va);

  return ret;
}
#endif

// -----------------------------------------------------------------------------
// Morphological operations for segmenting an image into text regions
// -----------------------------------------------------------------------------
static const char *segment_mask_sequence = "r11";
static const char *segment_seed_sequence = "r1143 + o4.4 + x4"; /* maybe o6.6 */
static const char *segment_dilation_sequence = "d3.3";

// -----------------------------------------------------------------------------
// Takes two pix as input, generated from the same original image:
//   1. pixb   - a binary thresholded image
//   2. piximg - a color or grayscale image
// and segments them by finding the areas that contain color or grayscale
// graphics.  These graphics regions are removed from the input binary
// image, and they are retained in the returned color-or-grayscale image.
// The upshot is that after this routine has been run:
//  (a) the input binary image contains only text, and is NULL if there
//      is no text, and
//  (b) the returned color-or-grayscale image contains only the graphics,
//      and is NULL if there is no graphics.
// The input color-or-grayscale image is not affected.
//
// Thanks to Dan Bloomberg for this
// -----------------------------------------------------------------------------

static PIX*
segment_image(PIX **ppixb, PIX *piximg)
{
  PIX *pixb = *ppixb;
  // Make a mask over the non-text (graphics) part of the input 1 bpp image
  // Do this by making a seed and mask, and filling the seed into the mask
  PIX *pixmask4 = pixMorphSequence(pixb, (char *) segment_mask_sequence, 0);
  PIX *pixseed4 = pixMorphSequence(pixb, (char *) segment_seed_sequence, 0);
  PIX *pixsf4 = pixSeedfillBinary(NULL, pixseed4, pixmask4, 8);
  PIX *pixd4 = pixMorphSequence(pixsf4, (char *) segment_dilation_sequence, 0);
  PIX *pixd = pixExpandBinaryPower2(pixd4, 4);
  pixDestroy(&pixd4);
  pixDestroy(&pixsf4);
  pixDestroy(&pixseed4);
  pixDestroy(&pixmask4);
  if (verbose) pixInfo(pixd, "mask image: ");

  // Remove pixels over the graphics part from the text mask.  This
  // side-effects the input binary mask.
  pixSubtract(pixb, pixb, pixd);

  // Set up table to count pixels in the text and graphics masks
  static l_int32 *tab = NULL;
  if (tab == NULL) tab = makePixelSumTab8();

  // If no graphics portion is found, destroy the graphics mask and return NULL
  l_int32  pcount;
  pixCountPixels(pixd, &pcount, tab);
  if (verbose) fprintf(stderr, "pixel count of graphics image: %u\n", pcount);
  if (pcount < 100)
  {
    pixDestroy(&pixd);
    return NULL;
  }

  // If no text portion is found, destroy the input binary image.
  pixCountPixels(pixb, &pcount, tab);
  if (verbose) fprintf(stderr, "pixel count of binary image: %u\n", pcount);
  if (pcount < 100)
  {
    pixDestroy(ppixb);  // destroy & set caller handle to NULL
    pixb = NULL;  // needed later in this function for pixInfo()
  }

  PIX *piximg1;
  if (piximg->d == 1 || piximg->d == 8 || piximg->d == 32)
  {
    piximg1 = pixClone(piximg);
  }
  else if (piximg->d > 8)
  {
    piximg1 = pixConvertTo32(piximg);
  }
  else
  {
    piximg1 = pixConvertTo8(piximg, FALSE);
  }

  PIX *pixd1;
  if (piximg1->d == 32)
  {
    pixd1 = pixConvertTo32(pixd);
  }
  else if (piximg1->d == 8)
  {
    pixd1 = pixConvertTo8(pixd, FALSE);
  }
  else
  {
    pixd1 = pixClone(pixd);
  }
  pixDestroy(&pixd);

  if (verbose)
  {
    pixInfo(pixd1, "binary mask image:");
    pixInfo(piximg1, "graphics image:");
  }
  pixRasteropFullImage(pixd1, piximg1, PIX_SRC | PIX_DST);

  pixDestroy(&piximg1);
  if (verbose)
  {
    pixInfo(pixb, "segmented binary text image:");
    pixInfo(pixd1, "segmented graphics image:");
  }

  return pixd1;
}

static int
get_ext_delim_pos(const char *fname)
{
  unsigned int pos = strcspn(fname,".");
  unsigned int last = 0;

  while (last + pos != strlen(fname))
  {
    last += (pos + 1);
    pos = strcspn(fname + last,".");
  }
  return last;
}

static char*
replace_suffix(char *name, const char *suffix)
{
  int extpos = get_ext_delim_pos(name);
  char *ret;
    
  asprintf(&ret, "%s ", name);
    
  ret[extpos] = '\0';;
  strcat(ret, suffix);
  return ret;
}

static char*
get_page_or_dict_name(char **elements, int cnt, const char *fname, int m_tiff)
{
  int i, extpos, same=-1;
  char *page_name, *pattern;
  const char *page_ext = ".jbig2";

  extpos = get_ext_delim_pos(fname);
  page_name = (char *) malloc(extpos + 12);
  memset(page_name,'\0',extpos + 12);
  if (extpos > 0)
    strncpy(page_name, fname, extpos-1);
  
  // Make sure first page from a multipage tiff also has a numerical extension
  // (this is necessary to guarantee it is sorted first)
  if (m_tiff)
    strcat(page_name, "_0000");
  strcat(page_name, page_ext);

  for (i=0; i<cnt; i++ )
  {
    if (strcmp(page_name, elements[i]) == 0)
    {
        same = i;
        break;
    }
  }

  if (same != -1)
  {
    int previdx=0, idx=0, res;

    pattern = (char *) malloc(extpos + 12);
    strcpy(pattern, page_name);
    pattern[extpos-1] = '\0';
    strcat(pattern, "_%4d.");

    for (i=same; i<cnt; i++)
    {
      res = sscanf(elements[i],pattern,&idx);
      if (res && idx > previdx) previdx = idx;
    }
    if (idx == 9999)
    {
      fprintf(stderr, "Cannot generate a unique name for %s\n", fname);
      exit(1);
    }
    sprintf(page_name + (extpos - 1), "_%04d%s", idx+1, page_ext);
    free(pattern);
    fprintf(stderr,"name %s\n",page_name);
  }
  return(page_name);
}

static int
is_tiff_format(int filetype)
{
  return (filetype == IFF_TIFF || 
          filetype == IFF_TIFF_PACKBITS || filetype == IFF_TIFF_RLE ||
          filetype == IFF_TIFF_G3 || filetype == IFF_TIFF_G4 ||
          filetype == IFF_TIFF_LZW || filetype == IFF_TIFF_ZIP);
}

int
main(int argc, char **argv)
{
  bool duplicate_line_removal = false;
  bool pdfmode = false;
  bool globalmode = false;
  int bw_threshold = BW_LOCAL_THRESHOLD_DEF;
  float threshold = JBIG2_THRESHOLD_DEF;
  float weight = JBIG2_WEIGHT_DEF;
  int pages_per_dict = 10;
  bool symbol_mode = false;
  bool refine = false;
  bool up2 = false, up4 = false;
  const char *output_threshold_image = NULL;
  l_int32 img_fmt = IFF_PNG;
  const char *img_ext = "bg.png";
  char **fnames;
  bool segment = false;
  bool auto_thresh = false;
  bool hash = true;
  int dpi = 0;
  int i, j;

  #ifdef WIN32
    int result = _setmode(_fileno(stdout), _O_BINARY);
    if (result == -1)
      fprintf(stderr, "Cannot set mode to binary for stdout\n");
  #endif

  for (i = 1; i < argc; ++i)
  {
    if (strcmp(argv[i], "-h") == 0 ||
        strcmp(argv[i], "--help") == 0)
    {
      usage(argv[0]);
      return 0;
      continue;
    }

    if (strcmp(argv[i], "-V") == 0 ||
        strcmp(argv[i], "--version") == 0)
    {
      fprintf(stderr, "jbig2enc %s\n", getVersion());
      return 0;
    }

    if (strcmp(argv[i], "-d") == 0 ||
        strcmp(argv[i], "--duplicate-line-removal") == 0)
    {
      duplicate_line_removal = true;
      continue;
    }

    if (strcmp(argv[i], "-p") == 0 ||
        strcmp(argv[i], "--pdf") == 0)
    {
      pdfmode = true;
      continue;
    }

    if (strcmp(argv[i], "-s") == 0 ||
        strcmp(argv[i], "--symbol-mode") == 0)
    {
      symbol_mode = true;
      continue;
    }

    if (strcmp(argv[i], "-r") == 0 ||
        strcmp(argv[i], "--refine") == 0)
    {
      fprintf(stderr, "Refinement broke in recent releases since it's "
                      "rarely used. If you need it you should bug "
                      "agl@imperialviolet.org to fix it\n");
      return 1;
      refine = true;
      continue;
    }

    if (strcmp(argv[i], "-2") == 0)
    {
      up2 = true;
      continue;
    }
    if (strcmp(argv[i], "-4") == 0)
    {
      up4 = true;
      continue;
    }

    if (strcmp(argv[i], "-O") == 0)
    {
      output_threshold_image = argv[i+1];
      i++;
      continue;
    }

    if (strcmp(argv[i], "-S") == 0)
    {
      segment = true;
      continue;
    }

    if (strcmp(argv[i], "-j") == 0 ||
        strcmp(argv[i], "--jpeg-output") == 0)
    {
      img_ext = "bg.jpg";
      img_fmt = IFF_JFIF_JPEG;
      continue;
    }

    if (strcmp(argv[i], "-t") == 0)
    {
      char *endptr;
      threshold = strtod(argv[i+1], &endptr);
      if (*endptr)
      {
        fprintf(stderr, "Cannot parse float value: %s\n", argv[i+1]);
        usage(argv[0]);
        return 1;
      }

      if ((threshold < JBIG2_THRESHOLD_MIN) ||
          (threshold > JBIG2_THRESHOLD_MAX))
      {
        fprintf(stderr, "Invalid value for threshold\n");
        fprintf(stderr, "(must be between %0.2f and %0.2f)\n",
                JBIG2_THRESHOLD_MIN, JBIG2_THRESHOLD_MAX);
        return 10;
      }
      i++;
      continue;
     }

    if (strcmp(argv[i], "-w") == 0)
    {
      char *endptr;
      weight = strtod(argv[i+1], &endptr);
      if (*endptr)
      {
        fprintf(stderr, "Cannot parse float value: %s\n", argv[i+1]);
        usage(argv[0]);
        return 1;
      }

      if ((weight < JBIG2_WEIGHT_MIN) || (weight > JBIG2_WEIGHT_MAX))
      {
        fprintf(stderr, "Invalid value for weight\n");
        fprintf(stderr, "(must be between %0.2f and %0.2f)\n",
                JBIG2_WEIGHT_MIN, JBIG2_WEIGHT_MAX);
        return 10;
      }
      i++;
      continue;
    }

    if (strcmp(argv[i], "-P") == 0 ||
        strcmp(argv[i], "--pages-per-dict") == 0)
    {
      char *endptr;
      pages_per_dict = strtol(argv[i+1], &endptr, 10);
      if (*endptr)
      {
        fprintf(stderr, "Cannot parse int value: %s\n", argv[i+1]);
        usage(argv[0]);
        return 1;
      }
      i++;
      continue;
    }

    // Local BW thresholding is the default.  However, if global
    // BW thresholding is requested, use its default threshold.
    if (strcmp(argv[i], "-G") == 0 ||
        strcmp(argv[i], "--global") == 0)
    {
      globalmode = true;
      bw_threshold = BW_GLOBAL_THRESHOLD_DEF;
      continue;
    }

    // If a BW threshold value is requested, overwrite the default value.
    if (strcmp(argv[i], "-T") == 0)
    {
      char *endptr;
      bw_threshold = strtol(argv[i+1], &endptr, 10);
      if (*endptr)
      {
        fprintf(stderr, "Cannot parse int value: %s\n", argv[i+1]);
        usage(argv[0]);
        return 1;
      }
      if (bw_threshold < BW_THRESHOLD_MIN || bw_threshold > BW_THRESHOLD_MAX)
      {
        fprintf(stderr, "Invalid bw threshold: (%d..%d)\n",
                BW_THRESHOLD_MIN, BW_THRESHOLD_MAX);
        return 11;
      }
      i++;
      continue;
    }

    // engage auto thresholding
    if (strcmp(argv[i], "--auto-thresh") == 0 ||
        strcmp(argv[i], "-a") == 0 )
    {
      auto_thresh = true;
      continue;
    }

    if (strcmp(argv[i], "--no-hash") == 0)
    {
      hash = false;
      continue;
    }

    if (strcmp(argv[i], "-v") == 0)
    {
      verbose = true;
      continue;
    }

    if (strcmp(argv[i], "-D") == 0 ||
        strcmp(argv[i], "--dpi") == 0)
    {
      char *endptr;
      long t_dpi = strtol(argv[i+1], &endptr, 10);
      if (*endptr)
      {
        fprintf(stderr, "Cannot parse int value: %s\n", argv[i+1]);
        usage(argv[0]);
        return 1;
      }
      if (t_dpi <= 0 || t_dpi > 9600)
      {
        fprintf(stderr, "Invalid dpi: (1..9600)\n");
        return 12;
      } 
      dpi = (int)t_dpi;
      i++;
      continue;
    }

    break;
  }

  if (i == argc)
  {
    fprintf(stderr, "No filename given\n\n");
    usage(argv[0]);
    return 4;
  }

  if (refine && !symbol_mode)
  {
    fprintf(stderr, "Refinement makes not sense unless in symbol mode!\n");
    fprintf(stderr, "(if you have -r, you must have -s)\n");
    return 5;
  }

  if (up2 && up4)
  {
    fprintf(stderr, "Can't have both -2 and -4!\n");
    return 6;
  }

  int numsubimages, subimage, num_pages = 0, num_images = argc - i, cnt = 0;
  char **images = argv + i;
  for (i = 0; i < num_images; i++)
  {
    subimage = numsubimages = 0;
    FILE *fp;
    if ((fp=lept_fopen(images[i], "r"))==NULL)
    {
      fprintf(stderr, "Unable to open \"%s\"\n", images[i]);
      return 1;
    }

    l_int32 filetype;
    findFileFormatStream(fp, &filetype);
    if (is_tiff_format(filetype) && tiffGetCount(fp, &numsubimages))
    {
      return 1;
    }
    fclose(fp);
    
    if (numsubimages) num_pages += numsubimages;
    else num_pages++;
  }
  subimage = numsubimages = i = 0;
  if (pages_per_dict <= 0) pages_per_dict = num_pages;
  fnames = (char **) malloc(sizeof(char *) * num_pages);

  for (i = cnt = 0; i < num_images && cnt < num_pages; )
  {
    struct jbig2ctx *ctx = jbig2_init(threshold, weight, 0, 0,
                                      !pdfmode, refine ? 10 : -1);
    int pages_to_compress = num_pages - cnt;
    int first_in_group = cnt;
    
    if (pages_to_compress > pages_per_dict)
      pages_to_compress = pages_per_dict;
    
    for (j = 0; j < pages_to_compress; cnt++, j++ )
    {
      if (subimage==numsubimages)
      {
        subimage = numsubimages = 0;
        FILE *fp;
        if ((fp=fopen(images[i], "r"))==NULL)
        {
          fprintf(stderr, "Unable to open \"%s\"\n", images[i]);
          return 1;
        }
        l_int32 filetype;
        findFileFormatStream(fp, &filetype);
        if (is_tiff_format(filetype) && tiffGetCount(fp, &numsubimages))
        {
          return 1;
        }
        fclose(fp);
      }

      PIX *source;
      if (numsubimages==0)
      {
        source = pixRead(images[i]);
      }
      else
      {
        source = pixReadTiff(images[i], subimage++);
      }

      if (dpi != 0 && source->xres == 0 && source->yres == 0)
      {
        source->xres = dpi;
        source->yres = dpi;
      }
      fnames[cnt] = get_page_or_dict_name(fnames, cnt, images[i], numsubimages > 1);

      if (!source) return 3;
      if (verbose)
        pixInfo(source, "source image:");

      PIX *pixl, *gray, *adapt, *pixt;
      if ((pixl = pixRemoveColormap(source, REMOVE_CMAP_BASED_ON_SRC)) == NULL)
      {
        fprintf(stderr, "Failed to remove colormap from %s\n", argv[i]);
        return 1;
      }
      pixDestroy(&source);

      if (pixl->d > 1)
      {
        if (pixl->d > 8)
        {
          gray = pixConvertRGBToGrayFast(pixl);
          if (!gray) return 1;
        }
        else if (pixl->d == 4 || pixl->d == 8)
        {
          gray = pixClone(pixl);
        }
        else
        {
          fprintf(stderr, "Unsupported input image depth: %d\n", pixl->d);
          return 1;
        }
        if (globalmode)
        {
          adapt = pixClone(gray);
        }
        else
        {
          adapt = pixCleanBackgroundToWhite(gray, NULL, NULL, 1.0, 90, 190);
        }
        pixDestroy(&gray);
        if (up2)
        {
          pixt = pixScaleGray2xLIThresh(adapt, bw_threshold);
        }
        else if (up4)
        {
          pixt = pixScaleGray4xLIThresh(adapt, bw_threshold);
        }
        else
        {
          pixt = pixThresholdToBinary(adapt, bw_threshold);
        }
        pixDestroy(&adapt);
      }
      else
      {
        pixt = pixClone(pixl);
      }

      if (verbose)
        pixInfo(pixt, "thresholded image:");

      if (output_threshold_image)
      {
        pixWrite(output_threshold_image, pixt, IFF_PNG);
      }

      if (segment && pixl->d > 1)
      {
        PIX *graphics = segment_image(&pixt, pixl);
        if (graphics)
        {
          if (verbose)
            pixInfo(graphics, "graphics image:");
          char *filename = replace_suffix(fnames[cnt], img_ext);
          pixWrite(filename, graphics, img_fmt);
          free(filename);
          pixDestroy(&graphics);
        }
        else if (verbose)
        {
          fprintf(stderr, "%s: no graphics found in input image\n", argv[i]);
        }
        if (!pixt)
        {
          fprintf(stderr, "%s: no text portion found in input image\n", argv[i]);
          i++;
          continue;
        }
      }
      pixDestroy(&pixl);

      if (!symbol_mode)
      {
        int length;
        uint8_t *ret;
        ret = jbig2_encode_generic(pixt, !pdfmode, 0, 0, duplicate_line_removal,
                                   &length);
        write(1, ret, length);
        return 0;
      }

      jbig2_add_page(ctx, pixt);
      pixDestroy(&pixt);
      if (subimage==numsubimages)
      {
        i++;
      }
    }

    if (auto_thresh)
    {
      if (hash)
      {
        jbig2enc_auto_threshold_using_hash(ctx);
      }
      else
      {
        jbig2enc_auto_threshold(ctx);
      }
    }

    uint8_t *ret;
    int length;
    ret = jbig2_pages_complete(ctx, &length);
    if (pdfmode)
    {
      char *dict_name = replace_suffix(fnames[first_in_group], "sym");
      const int fd = open(dict_name, O_WRONLY | O_TRUNC | O_CREAT | WINBINARY, 0600);
      free(dict_name);
      if (fd < 0) abort();
      write(fd, ret, length);
      close(fd);
    }
    else
    {
      write(1, ret, length);
    }
    free(ret);

    for (j = 0; j < pages_to_compress; ++j)
    {
      ret = jbig2_produce_page(ctx, j, -1, -1, &length);
      if (pdfmode)
      {
        const int fd = open(fnames[cnt - pages_to_compress + j], O_WRONLY | O_CREAT | O_TRUNC | WINBINARY, 0600);
        if (fd < 0) abort();
        write(fd, ret, length);
        close(fd);
      }
      else
      {
        write(1, ret, length);
      }
      free(ret);
    }
    jbig2_destroy(ctx);
  }
  for (i=0; i<cnt; i++) free(fnames[i]);
  free(fnames);

  return 0;
}
