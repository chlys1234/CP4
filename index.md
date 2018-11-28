# Computed Photography Assignment 4

### HDR Imaging
-------------
I was able to import a given .NEF format raw file into a 16bit-TIFF file on MATLAB using dcraw. As given in the assignment,
- Do white balancing using the camera’s profile for white balancing
- Do demosaicing using high-quality interpolation
- Use sRGB as the output color space

After reading the dcraw documentation, I loaded a .NEF format file from MATLAB with the following flag:
``` im = imread(dc,'exposure10.nef','-w -T -6 -q 3');```

### Linearize Rendered Images
-------------
Unlike a raw file, a rendered image has a non-linear property, so it needs to be linearized. To do this, we refer to Debevec’s paper Recovering high dynamic range radiance maps from photographs and his code. To perform the linearize, he assumed that the pixel of the image is associated with unknown space radiance value and shutter speed. 
