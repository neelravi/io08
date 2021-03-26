#!/bin/bash
FC=ifort

$FC -c m_periodic_table.F90 m_keywords.F90
$FC interface.F90 m_keywords.o m_periodic_table.o /usr/local/lib/libfdf.a
