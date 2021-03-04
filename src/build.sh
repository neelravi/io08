#!/bin/bash
ifort -c periodic_table_m.F90
ifort -c keywords_m.F90

ifort iochamp.f90 keywords_m.o periodic_table_m.o /usr/local/lib/libfdf.a
