# A Modern Fortran-based Parser

<!-- Thanks for visiting [The Markdown Guide](https://www.markdownguide.org)!

This Markdown cheat sheet provides a quick overview of all the Markdown syntax elements. It can’t cover every edge case, so if you need more information about any of these elements, refer to the reference guides for [basic syntax](https://www.markdownguide.org/basic-syntax) and [extended syntax](https://www.markdownguide.org/extended-syntax). -->

## Get the code
  The parser uses a modified libfdf library. This is included in this repository as a submodule. To clone the entire project, do

  `git clone --recurse-submodules https://github.com/TREX-CoE/iof08.git`


## Compilation
  The project contains two directories (a) modified-libfdf and (b) parser.

  - Compile and install the modified-libfdf using the following set of commands
    - `./configure --prefix=/usr/local FC=ifort CC=icc `
    - `make`
    - `sudo make install`
    - `make check`


## Integrate parser in your code



### Features of the parser (including inheritance from libfdf)

- Include another input file for parser to read using:

` %include    global.inp`

- Include a data file for parser to read using:

` load label filename`

Here, depending upon the label, parser will provide the filename. For example,

` load basis cc-pvtz.gbs`

- Read molecular coordinates directly from the input file using 

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

- Read molecular coordinates from an external .xyz file using 

` %block molecule < caffeine.xyz `


- Group certain keywords using the %module construct

```perl
%module DMC
    tau     =   0.04
    etrial  = -15 Ha
%endmodule
```

- Logical variables accept true, TRUE, T, 1, .true. as valid keywords for true
` optimize_wavefunction 	true`




**bold text**

### Italic

*italicized text*

### Blockquote

> blockquote

### Ordered List

1. First item
2. Second item
3. Third item

### Unordered List

- First item
- Second item
- Third item

### Code

`code`

### Horizontal Rule

---

### Link

[title](https://www.example.com)

### Image

![alt text](image.jpg)

## Extended Syntax

These elements extend the basic syntax by adding additional features. Not all Markdown applications support these elements.

### Table

| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |

### Fenced Code Block

```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```

### Footnote

Here's a sentence with a footnote. [^1]

[^1]: This is the footnote.

### Heading ID

### My Great Heading {#custom-id}

### Definition List

term
: definition

### Strikethrough

~~The world is flat.~~

### Task List

- [x] Write the press release
- [ ] Update the website
- [ ] Contact the media
