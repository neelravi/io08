#!/bin/bash
FC=ifort
echo $FC

$FC -c m_periodic_table.F90
$FC -c m_keywords.F90

$FC iochamp.f90 m_keywords.o m_periodic_table.o /usr/local/lib/libfdf.a
