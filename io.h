#ifndef IO_H
#define IO_H
//#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
//ldoc on
/**
 * ## Particle I/O
 *
 * Writing out the particle positions is expensive, but it seems
 * necessary if we want to make pretty pictures of how they are all
 * moving.  The particle I/O system writes out a simple text-based CSV
 * (comma-separated value) file with the fields:
 *
 * - `PTag`: Indicates whether particles is black (0) or red (1)
 * - `PId`: Unique integer identifier for the particle
 * - `PLocX`, `PLocY`: Position components
 * - `PDirX`, `PDirY`: Velocity components
 *
 * The CSV file may have several frames; they are assumed to be in
 * consecutive order.
 */

FILE* start_frames(const char* fname, int npart);
void write_frame(FILE* fp, float* x, float* rad, int* type, int n, float L, float* col);
void end_frames(FILE* fp);

#endif
//ldoc off
