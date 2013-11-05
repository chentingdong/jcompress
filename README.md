       
                   Project for course CS175 (Dater Compressing)
              
                                            By Tingdong Chen
                                            Dec 4, 2000.

[Joking]
CTD compresser and decompresser for ppm format images.
VERSION 0.01.

[Usage]

1. To compile:

>g++ -o ctd ctd.cpp

2. To use:

>ctd test.ppm

It will automatically give out a compressed file called encode.ctd, and a decompressed 
file from encode.ctd called decode.ppm.

3. To view images:

I did all this project under Linux Red Hat 6.1, where I can view the ppm files with
xv. Don't know how to view it under windows.

4. change compression rate (and thus the quality)

change the var 'quality' in the source codes, change it larger will make the decode.ppm 
looks better, but the encode.ctd larger.


[Comparation]

1. The roses in color I hand in is a comparation between my decode.ppm with the 
original test.ppm. It's enlarged nearly 10 times of the original size just to let you 
see the details. The gray pair are the original size.

2. size comparation: (with approximately the same quality)

 test.ppm      412656
 encode.ctd     30658
 test.jpg        7962

3. time consumed:
 not very fast, this rose I used nearly 2 seconds on my PIII500.

[Known bugs and later improvement:]

1. for some image, the white part has some minus value, then the xv can't display 
the decode.ppm. I simply changed the minus value to 0; which is clearly seen in my
sample of decode.ppm. I didn't have time to correct it, but I think it's not hard.

2. for ppm files that has notes begin with '#', I forgot to tell the program to 
skip it, just delete that line.

3. for very small images, after decompression, I can clearly see the sudden change
between neighbor blocks. I should do some smoothing between neighbor 8X8 blocks later.

4. It's very difficult for me to understand the format of images, that's why I choose
the ppm format, which is simply ascii. I will, if free, extend it for some other format
images, such as bmp gif....

Thank you. 
