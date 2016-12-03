//
//  io_file.h
//  demo_io_file
//
//  Created by KaiWu on Aug/10/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef io_file_h
#define io_file_h

#include <stdio.h>
#include <stdlib.h>
#include "RPly/rply.h"
#include "RPly/rplyfile.h"

#define DIM 3

p_ply ply_header_read(const char *fname, int *num_vertices, int *num_faces, int *num_range_grid);

int ply_read_1(const char *fname, double *pt, double *n, double *rgb, double *face, double *grid);

int ply_write_1(const char *fname, double *pt, double *n, double *rgb, double *face, double *grid, int num_vertices, int num_faces, int row, int col);

int ply_write_list(const char *fname, double *grid, int row, int col);

int txt_read();

int txt_write();

int obj_read(const char *fname, double *pt);

int obj_write(const char *fname, double *pt);

static int read_scalar_cb(p_ply_argument argument);

static int read_list_cb(p_ply_argument argument);

// helper methods
/**
 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C":
 */

void strreverse(char* begin, char* end);
void myitoa(int value, char* str, int base);


#endif /* io_file_h */
