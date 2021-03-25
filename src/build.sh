#!/bin/bash
ifort -c m_periodic_table.F90
ifort -c m_keywords.F90

ifort iochamp.f90 m_keywords.o m_periodic_table.o /usr/local/lib/libfdf.a
