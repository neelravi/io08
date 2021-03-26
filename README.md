# A Modern Fortran-based Parser

  This repository contains a modern Fortran-based input file parser. It uses a modified libfdf library.

---

## Get the code
  The library is included in this repository as a submodule. To clone the entire project, do

  `git clone --recurse-submodules https://github.com/TREX-CoE/fparser.git`


## Compilation
  The project contains two folders (a) modified-libfdf and (b) parser.

  Compile and install the modified-libfdf using the following set of commands

  `./configure --prefix=/usr/local FC=ifort CC=icc`

  `make`

  `sudo make install`

  In the parser folder, link the modified libfdf library with the interface Fortran file.

  `ifort -c interface.F90 m_keywords.F90 m_periodic_table.F90 /usr/local/lib/libfdf.a`


## Integrate parser in your code
  Just include the interface.F90 file and the keyword declaration module files in the existing Makefile of your code.

---

## Features of the parser (including inheritance from libfdf)

- Blank lines and lines starting with # are ignored.

- The order of keyword-value pairs does not matter

- Multiple keyword-value pairs can be put in a single line

- Keywords are case insensitive.

- A default value can be set for keywords not present in the input file. 


## Syntax

1. Include another input file for parser to read using:

    ` %include    global.inp`

2. Include a data file for parser to read using:

    ` load label filename`

3. Here, depending upon the label, parser will provide the filename. For example,

    ` load basis cc-pvtz.gbs`

4. Read molecular coordinates directly from the input file using 

    ```perl
    %block molecule 
    12
    #benzene comment
    C    0.00000    1.40272  0
    H    0.00000    2.49029  0
    C   -1.21479    0.70136  0
    H   -2.15666    1.24515  0
    C   -1.21479   -0.70136  0
    H   -2.15666   -1.24515  0
    C    0.00000   -1.40272  0
    H    0.00000   -2.49029  0
    C    1.21479   -0.70136  0
    H    2.15666   -1.24515  0
    C    1.21479    0.70136  0
    H    2.15666    1.24515  0
    %endblock
    ```

5. Read molecular coordinates from an external .xyz file using 

    ` %block molecule < benzene.xyz `

6. Group certain keywords using the %module construct

    ```perl
    %module DMC
      tau     =   0.04
      etrial  = -15 Ha
    %endmodule
    ```

7. Logical variables accept `true`, `TRUE`, `T`, `.true.` as valid keywords for `True`. The `fdf_boolean` function can also  take "1" as true and "0" as false from the input.

    ` optimize_wavefunction 	true`

8. Single and Double precision numbers along with numbers in scientific format can be read using the `fdf_get()` function.

    `energy_tol = fdf_get('energy_tol', 0.00001d0)`

9. Floats/integers/strings/booleans can be parsed generically using the interface `fdf_get()` function. Strings are limited to 132 characters per line.

    `  sr_tau         = fdf_get('sr_tau', 0.025d0)`

    `  nspin1         = fdf_get('nspin1', 1)`

    `  opt_method     = fdf_get('opt_method', "sr_n")`

    `  multiple_adiag = fdf_get('multiple_adiag', .false.)`

10. Units can be specified to variables. Unit conversion is possible at the parsing.

    If the input file has `etrial  = -15 Ha` entry, the `etrial` variable can be assigned values using the `fdf_physical` function with unit conversion.

    `  etrial = fdf_physical('etrial', -20.d0, 'eV')`

11. List of public functions available for parsing the data:

    Initiate

    `fdf_init` 

    `fdf_shutdown`

    Single data

    `fdf_get`

    `fdf_integer`

    `fdf_single`

    `fdf_double`

    `fdf_string`

    `fdf_boolean`

    `fdf_physical`

    `fdf_convfac`

    `fdf_load_filename`

    Lists (data enclosed in [])

    `fdf_islist`

    `fdf_islinteger`

    `fdf_islreal`

    `fdf_list`

    `fdf_linteger`

    `fdf_ldouble`

    Returns the string associated with a mark line

    `fdf_getline`

    Test if a label is defined

    `fdf_defined`

    `fdf_isphysical`

    `fdf_isblock`

    `fdf_load_defined`

    Allow to overwrite things in the FDF

    `fdf_overwrite`

    `fdf_removelabel`

    `fdf_addline`

    Block reading (processing each line of data)

    `fdf_block, fdf_block_linecount`

    `fdf_bline, fdf_bbackspace, fdf_brewind, fdf_bclose`

    `fdf_bnintegers, fdf_bnreals, fdf_bnvalues, fdf_bnnames, fdf_bntokens`

    `fdf_bintegers, fdf_breals, fdf_bvalues, fdf_bnames, fdf_btokens`

    `fdf_bboolean, fdf_bphysical`

    `fdf_bnlists, fdf_bnilists, fdf_bnrlists, fdf_bnvlists`

    `fdf_bilists, fdf_brlists, fdf_bvlists`
  
    Match, search over blocks, and destroy block structure

    `fdf_bmatch, fdf_bsearch, fdf_substring_search`

    `fdf_setoutput, fdf_setdebug`

---

## Demonstration


  In the `parser` folder, we have included a sample `interface.F90` and `m_keywords.F90` files. 
  In the `interface.F90` file, we have demonstrated how keyword-values pairs, simple data blocks, 
  and data from external files can be read easily. 

