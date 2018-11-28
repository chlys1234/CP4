# Computed Photography Assignment 4

### HDR Imaging
I was able to import a given .NEF format raw file into a 16bit-TIFF file on MATLAB using dcraw. As given in the assignment,
- Do white balancing using the camera’s profile for white balancing
- Do demosaicing using high-quality interpolation
- Use sRGB as the output color space
After reading the dcraw documentation, I loaded a .NEF format file from MATLAB with the following flag:
