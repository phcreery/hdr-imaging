#pragma once
#include <stdint.h>

/*
 * Load PNG/JPG images from explicit paths, parse exposure times from filenames,
 * sort by exposure descending, and return the log-exposure array B.
 *
 * Filename format: {numerator}_{denominator}.<ext>
 *   e.g. "32_1.png"  -> exposure = 32.0 s
 *        "1_4.png"   -> exposure = 0.25 s
 *
 * Returns:
 *   Heap-allocated array of num_paths pointers, each pointing to a
 *   W*H*3 uint8_t buffer (RGB, row-major interleaved).
 *   B_out[i] = log(exposure) for sorted image i (descending exposure).
 *   *w_out and *h_out are set to image dimensions (all images must match).
 *
 * Caller must free each pixel buffer and the pointer array via free_images().
 * Returns NULL on failure.
 */
uint8_t **load_images(const char **paths, int num_paths,
                      int *w_out, int *h_out,
                      float *B_out);

void free_images(uint8_t **images, int num_paths);
