//
//  io_file.c
//  demo_io_file
//
//  Created by KaiWu on Aug/10/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "io_file.h"

double *point;
double *normal;
double *rgb;
double *face;
double *grid;
int type_property[5]; // 1: pt, 2: n, 3: rgb, 4: face, 5: range_grid
int ind_property = 0;
int num_property = 1; // 1: pt, 2: pt and n, or pt and rgb, 3: pt, n, and rgb

p_ply ply_header_read(const char *fname, int *num_vertices, int *num_faces, int *num_range_grid)
{
    p_ply ply = ply_open(fname, NULL, 0, NULL);
    
    if(!ply)
        return NULL;
    
    if(!ply_read_header(ply))
        return NULL;
    
    p_ply_element element = NULL;
    while ((element = ply_get_next_element(ply, element)))
    {
        p_ply_property property = NULL;
        long ninstances;
        const char *element_name;
        ply_get_element_info(element, &element_name, &ninstances);
        if(!strcmp(element_name, "vertex"))
        {
            *num_vertices = (int)ninstances;
            
            while((property = ply_get_next_property(element, property)))
            {
                const char *property_name;
                ply_get_property_info(property, &property_name, NULL, NULL, NULL);
                if(!strcmp(property_name, "x") || !strcmp(property_name, "y") || !strcmp(property_name, "z"))
                    type_property[0] = 1;
                else if(!strcmp(property_name, "nx") || !strcmp(property_name, "ny") || !strcmp(property_name, "nz"))
                    type_property[1] = 1;
                else if(!strcmp(property_name, "red") || !strcmp(property_name, "green") || !strcmp(property_name, "blue"))
                    type_property[2] = 1;
            }
        }
        else if(!strcmp(element_name, "face") && num_faces != NULL)
        {
            *num_faces = (int)ninstances;
            type_property[3] = 1;
        }
        else if(!strcmp(element_name, "range_grid") && num_range_grid != NULL)
        {
            *num_range_grid = (int)ninstances;
            type_property[4] = 1;
        }
    }
    
    return ply;
}

int ply_read_1(const char *fname, double *pt, double *n, double *c, double *f, double *g)
{
    int num_vertices;
    int num_grids;
    
    p_ply ply = ply_header_read(fname, &num_vertices, NULL, &num_grids);
    
    point = pt;
    normal = n;
    rgb = c;
    face = f;
    grid = g;
    
    ply_set_read_cb(ply, "vertex", "x", read_scalar_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "y", read_scalar_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "z", read_scalar_cb, NULL, 0);
    num_property = 1;
    
    if(n != NULL && c == NULL && type_property[1])
    {
        ply_set_read_cb(ply, "vertex", "nx", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "ny", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "nz", read_scalar_cb, NULL, 0);
        num_property = 2;
    }
    else if(n == NULL && c != NULL && type_property[2])
    {
        ply_set_read_cb(ply, "vertex", "red", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "green", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "blue", read_scalar_cb, NULL, 0);
        num_property = 2;
    }
    else if(n != NULL && c != NULL && type_property[1] && type_property[2])
    {
        ply_set_read_cb(ply, "vertex", "nx", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "ny", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "nz", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "red", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "green", read_scalar_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "blue", read_scalar_cb, NULL, 0);
        num_property = 3;
    }
    
    if(f != NULL && type_property[3])
    {
        ply_set_read_cb(ply, "face", "ind_vertex", read_list_cb, NULL, 0);
    }
    
    if(g != NULL && type_property[4])
    {
        ply_set_read_cb(ply, "range_grid", "vertex_indices", read_list_cb, NULL, 0);
    }
    
    if (!ply_read(ply))
        return -1;
    
    if(!ply_close(ply))
        return -1;
    
    return 0;
}

int read_scalar_cb(p_ply_argument argument)
{
    double val = ply_get_argument_value(argument);
    int ind_instance = (int)get_argument_instance_index(argument);
    
    switch(num_property)
    {
        case 1:
            point[DIM * ind_instance + ind_property] = val;
            
            ind_property++;
            ind_property = ind_property % DIM;
            break;
        case 2:
            if(ind_property < DIM)
                point[DIM * ind_instance + ind_property] = val;
            else if(normal != NULL && rgb == NULL)
                normal[DIM * ind_instance + ind_property - DIM] = val;
            else if(rgb != NULL)
                rgb[DIM * ind_instance + ind_property - DIM] = val;
            
            ind_property++;
            ind_property = ind_property % (2 * DIM);
            break;
        case 3:
            if(ind_property < DIM)
                point[DIM * ind_instance + ind_property] = val;
            else if(ind_property < 2 * DIM)
                normal[DIM * ind_instance + ind_property - DIM] = val;
            else
                rgb[DIM * ind_instance + ind_property - 2 * DIM] = val;
            
            ind_property++;
            ind_property = ind_property % (3 * DIM);
            break;
    }
    
    return 1; // restricted by the RPly library
}

int read_list_cb(p_ply_argument argument)
{
    long length, ind_val;
    ply_get_argument_property(argument, NULL, &length, &ind_val);
    
    int ind_instance = (int)get_argument_instance_index(argument);
    double val = ply_get_argument_value(argument);
    
    // For range_grid
    if(ind_val >= 0)
        grid[ind_instance] = val;
    else
        grid[ind_instance] = -1;
    
    // For face, no yet
    
    return 1; // restricted by the RPly library
}

int ply_write_list(const char *fname, double *g, int row, int col)
{
    p_ply ply = ply_create(fname, PLY_ASCII, NULL, 0, NULL);
    
    if (!ply)
        return -1;
    
    int num_grids = row * col;
    ply_add_element(ply, "range_grid", num_grids);
    ply_add_list_property(ply, "vertex_indices", PLY_UCHAR, PLY_INT);
    
    if(!ply_write_header(ply))
        return -1;
    
    for (size_t i = 0; i < num_grids; ++i)
    {
        if(g[i] >=0)
        {
            ply_write(ply, 1);
            ply_write(ply, g[i]);
        }
        else
            ply_write(ply, 0);
    }
    
    return 0;
}

// no callback function
int ply_write_1(const char *fname, double *pt, double *n, double *c, double *f, double *g, int num_vertices, int num_faces, int row, int col)
{
    p_ply ply = ply_create(fname, PLY_ASCII, NULL, 0, NULL);
    
    if(!ply)
        return -1;
    
    // add obj_info
    char *obj_info_col_0 = "num_cols ";
    char *obj_info_row_0 = "num_rows ";
    char *obj_info_col_1 = NULL, *obj_info_row_1 = NULL;
    if(col < 100)
        obj_info_col_1 = (char *)malloc(2 * 8);
    else if(col < 1000)
        obj_info_col_1 = (char *)malloc(3 * 8);
    else if(col < 10000)
        obj_info_col_1 = (char *)malloc(4 * 8);
    
    if(row < 100)
        obj_info_row_1 = (char *)malloc(2 * 8);
    else if(row  < 1000)
        obj_info_row_1 = (char *)malloc(3 * 8);
    else if(row < 10000)
        obj_info_row_1 = (char *)malloc(4 * 8);
    myitoa (col, obj_info_col_1, 10);
    myitoa (row, obj_info_row_1, 10);
    
    char *obj_info_col = (char *)malloc(strlen(obj_info_col_0) + strlen(obj_info_col_1) + 1);
    strcpy(obj_info_col, obj_info_col_0);
    strcat(obj_info_col, obj_info_col_1);
    char *obj_info_row = (char *)malloc(strlen(obj_info_row_0) + strlen(obj_info_row_1) + 1);
    strcpy(obj_info_row, obj_info_row_0);
    strcat(obj_info_row, obj_info_row_1);
    ply_add_obj_info(ply, obj_info_col);
    ply_add_obj_info(ply, obj_info_row);
    
    // add elements and properties
    ply_add_element(ply, "vertex", num_vertices);
    ply_add_scalar_property(ply, "x", PLY_FLOAT);
    ply_add_scalar_property(ply, "y", PLY_FLOAT);
    ply_add_scalar_property(ply, "z", PLY_FLOAT);
    if(n != NULL)
    {
        ply_add_scalar_property(ply, "nx", PLY_FLOAT);
        ply_add_scalar_property(ply, "ny", PLY_FLOAT);
        ply_add_scalar_property(ply, "nz", PLY_FLOAT);
    }
    if(c != NULL)
    {
        ply_add_scalar_property(ply, "red", PLY_UCHAR);
        ply_add_scalar_property(ply, "green", PLY_UCHAR);
        ply_add_scalar_property(ply, "blue", PLY_UCHAR);
    }
    if(f != NULL && num_faces > 0)
    {
        ply_add_element(ply, "face", num_faces);
        ply_add_list_property(ply, "vertex_indices", PLY_UCHAR, PLY_INT);
    }
    int num_grids = row * col;
    if(g != NULL && num_grids > 0)
    {
        ply_add_element(ply, "range_grid", num_grids);
        ply_add_list_property(ply, "vertex_indices", PLY_UCHAR, PLY_INT);
    }
    
    if (!ply_write_header(ply))
        return -1;
    
    // write vertex
    for (size_t i = 0; i < num_vertices; ++i)
    {
        ply_write(ply, pt[DIM * i + 0]);
        ply_write(ply, pt[DIM * i + 1]);
        ply_write(ply, pt[DIM * i + 2]);
        if (n != NULL)
        {
            ply_write(ply, n[DIM * i + 0]);
            ply_write(ply, n[DIM * i + 1]);
            ply_write(ply, n[DIM * i + 2]);
        }
        if (c != NULL)
        {
            ply_write(ply, c[DIM * i + 0]);
            ply_write(ply, c[DIM * i + 1]);
            ply_write(ply, c[DIM * i + 2]);
        }
    }
    
    // write faces, not checked
    if(f != NULL && num_faces > 0)
    {
        size_t nf = 0;
        size_t i = 0;
        while (nf < num_faces)
        {
            size_t nv = f[i];
            for (int j = 0; j <= nv; ++j) // write number of vertices and indices of vertices
            {
                ply_write(ply, f[i + j]);
            }
            i = i + nv;
            nf++;
        }
    }
    
    if(g != NULL && num_grids > 0)
    {
        for (size_t i = 0; i < num_grids; ++i)
        {
            if (g[i] >= 0)
            {
                ply_write(ply, 1);
                ply_write(ply, g[i]);
            }
            else
                ply_write(ply, 0);
        }
    }
    
    if(!ply_close(ply))
        return -1;
    
    return 0;
}

void strreverse(char* begin, char* end) {
    char aux;
    while(end>begin)
        aux=*end, *end--=*begin, *begin++=aux;
}

void myitoa(int value, char* str, int base) {
    static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char* wstr=str;
    int sign;
    
    // Validate base
    if (base<2 || base>35){ *wstr='\0'; return; }
    
    // Take care of sign
    if ((sign=value) < 0) value = -value;
    
    // Conversion. Number is reversed.
    do *wstr++ = num[value%base]; while(value/=base);
    if(sign<0) *wstr++='-';
    *wstr='\0';
    
    // Reverse string
    strreverse(str,wstr-1);
}
