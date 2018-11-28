# Computed Photography Assignment 4

### HDR Imaging
-------------
I was able to import a given .NEF format raw file into a 16bit-TIFF file on MATLAB using dcraw. As given in the assignment,
- Do white balancing using the camera’s profile for white balancing
- Do demosaicing using high-quality interpolation
- Use sRGB as the output color space

After reading the dcraw documentation, I loaded a .NEF format file from MATLAB with the following flag:

<pre><code> im = imread(dc,'exposure10.nef','-w -T -6 -q 3');
</code></pre>

### Linearize Rendered Images
-------------
Unlike a raw file, a rendered image has a non-linear property, so it needs to be linearized. To do this, we refer to Debevec’s paper Recovering high dynamic range radiance maps from photographs and his code. To perform the linearize, he assumed that the pixel of the image is associated with unknown space radiance value and shutter speed. 

![image](https://user-images.githubusercontent.com/45420635/49143846-46bbc480-f33f-11e8-84bc-61a7a4ec761f.png)

In order to make the pixel mapping as smooth as possible, the following optimization problem was solved through MATLAB with the g function. It is quite difficult because L term is unknown. But Debevec has kindly attached the gsolve function to solve it.

![image](https://user-images.githubusercontent.com/45420635/49143947-7f5b9e00-f33f-11e8-97ce-7b1bce5d78b8.png)

The code that Debevec put in the paper is as follows.

<pre><code> function [g,lE]=gsolve(Z,B,w)
n = 256;
A = sparse((size(Z,1)/100)*size(Z,2)+n+1,n+(size(Z,1)/100));
b = sparse(size(A,1),1);

%% Include the data fitting equations

k = 1;
for j=1:size(Z,2)
    for i=1:100:size(Z,1)
        wij = w(Z(i,j)+1);

        if Z(i,j) > 0
            A(k,Z(i,j)+1) = wij;
            A(k,n+i) = -wij;
            b(k,1) = wij * B(k,1);
        else
            A(k,Z(i,j)+1) = 0;
            b(k,1) = 0;
        end   
        k=k+1;
    end
end

%% Fix the curve by setting its middle value to 0

A(k,129) = 1;
k=k+1;

%% Include the smoothness equations

for i=1:n-2
    A(k,i)=1*w(i+1);
    A(k,i+1)=-2*w(i+1);
    A(k,i+2)=1*w(i+1);

    k=k+1;
end
%% Solve the system using SVD
x = A\b;
g = x(1:n);
lE = x(n+1:size(x,1));
</code></pre>

Since the image size was very large (4000 X 6000), we resized it 0.1x and in particular 100:1 sampling pixels into the equation. As a result, Pixel’s linear mapping looks like this. In the above equation, both the uniform type and the Tent type was used as a  w function as shown below.

<pre><code> w_uni = zeros(256,1);
w_uni(3:253) = 1./256;
w_ten = zeros(256,1);

for k=1:256
    if k <= 128
        w_ten(k) = (k)./256;
    else
        w_ten(k) = (257-k)./256;
    end
end
</code></pre>

![image](https://user-images.githubusercontent.com/45420635/49144071-c9448400-f33f-11e8-8a96-361052872847.png)

Here’s my code

<pre><code> M=400;
N=600;
K=16;

for i=1:16
    i
    filename = sprintf('exposure%d.jpg',i); 
    im = imread(filename);
    im = imresize(im,0.1);
    Z(:,:,:,i) = im; 
    B(1+(i-1)*(M/10)*(N/10):i*(M/10)*(N/10),1) = log((2^(i-1))/2048);
end

Z = reshape(Z,M*N,3,16);

for c=1:3
    [g,lE]=gsolve(squeeze(Z(:,c,:)),B,w_ten);
    g_ten(:,c) = g;
end
g_ten = full(g_ten);

for c=1:3
    [g,lE]=gsolve(squeeze(Z(:,c,:)),B,w_uni);
    g_uni(:,c) = g;
end
g_uni = full(g_uni);

for c=1:3
    c
    g = g_ten(:,c);
    Z_lin_ten(:,c,:) = exp(g(Z(:,c,:)+1)./255);
end

for k=1:16
    for c=1:3
        Z_lin_ten(:,c,k) = uint8(255*(Z_lin_ten(:,c,k)-min(min(min(Z_lin_ten(:,c,k)))))./(max(max(max(Z_lin_ten(:,c,k))))-min(min(min(Z_lin_ten(:,c,k))))));
    end
end
Z_temp = reshape(Z,M,N,3,16);

for c=1:3
    c
    g = g_uni(:,c);
    Z_lin_uni(:,c,:) = exp(g(Z(:,c,:)+1)./255);
end

for k=1:16
    for c=1:3
        Z_lin_uni(:,c,k) = uint8(255*(Z_lin_uni(:,c,k)-min(min(min(Z_lin_uni(:,c,k)))))./(max(max(max(Z_lin_uni(:,c,k))))-min(min(min(Z_lin_uni(:,c,k))))));
    end
end
</code></pre>

* Weighting function = Uniform

![image](https://user-images.githubusercontent.com/45420635/49144187-132d6a00-f340-11e8-9bf9-8bf9714d8a45.png)

* Weighting function = Tent

![image](https://user-images.githubusercontent.com/45420635/49144224-2dffde80-f340-11e8-97c0-7fde0a0ac81e.png)

Plot results show that the pixel mapping is smoother than uniform when the weighting function is tent.

### Merge Exposure stack into HDR Image
-------------
The next step is to pass the LDR images through the g function, then make each HDR image and Merge into one image. For JPEG image, we tried to perform linear merging and logarithm merging when the weighting function is uniform and when it is tent.
When using linear merging, the HDR image is formed as:

![image](https://user-images.githubusercontent.com/45420635/49144256-4cfe7080-f340-11e8-8737-de0a17d236d6.png)

When using logarithmic merging, the HDR image is formed as:

![image](https://user-images.githubusercontent.com/45420635/49144291-5e477d00-f340-11e8-9836-1065c1b6076e.png)

Here’s my code

<pre><code> Z_lin_ten = reshape(Z_lin_ten,M,N,3,16);
Z_lin_uni = reshape(Z_lin_uni,M,N,3,16);

% Linear merging

for k=1:16
    num_temp = w_ten(Z_temp(:,:,:,k)+1).*(Z_lin_ten(:,:,:,k)./255)./(2^(k-1)/2048);
    den_temp = w_ten(Z_temp(:,:,:,k)+1);
    num = num + num_temp;
    den = den + den_temp;
end

final_linm_tent = num./den;
for c=1:3
    final_linm_tent(:,:,c) = 255.*(final_linm_tent(:,:,c) - min(min(min(final_linm_tent(:,:,c)))))./(max(max(max(final_linm_tent(:,:,c))))-min(min(min(final_linm_tent(:,:,c)))));
end

for k=1:16
    num_temp = w_uni(Z_temp(:,:,:,k)+1).*(Z_lin_uni(:,:,:,k)./255)./(2^(k-1)/2048);
    den_temp = w_uni(Z_temp(:,:,:,k)+1);
    num = num + num_temp;
    den = den + den_temp;
end

final_linm_uni = num./den;
for c=1:3
    final_linm_uni(:,:,c) = 255.*(final_linm_uni(:,:,c) - min(min(min(final_linm_uni(:,:,c)))))./(max(max(max(final_linm_uni(:,:,c))))-min(min(min(final_linm_uni(:,:,c)))));
end

% Logarithmic merging

for k=1:16
    num_temp = w_ten(Z_temp(:,:,:,k)+1).*(Z_lin_ten(:,:,:,k)./255-log(2^(k-1)/2048));
    den_temp = w_ten(Z_temp(:,:,:,k)+1);
    num = num + num_temp;
    den = den + den_temp;
end

final_logm_tent = exp(num./den./255);
for c=1:3
    final_logm_tent(:,:,c) = 255.*(final_logm_tent(:,:,c) - min(min(min(final_logm_tent(:,:,c)))))./(max(max(max(final_logm_tent(:,:,c))))-min(min(min(final_logm_tent(:,:,c)))));
end

for k=1:16
    num_temp = w_uni(Z_temp(:,:,:,k)+1).*(Z_lin_uni(:,:,:,k)./255-log(2^(k-1)/2048));
    den_temp = w_uni(Z_temp(:,:,:,k)+1);
    num = num + num_temp;
    den = den + den_temp;
end

final_logm_uni = exp(num./den./255);
for c=1:3
    final_logm_uni(:,:,c) = 255.*(final_logm_uni(:,:,c) - min(min(min(final_logm_uni(:,:,c)))))./(max(max(max(final_logm_uni(:,:,c))))-min(min(min(final_logm_uni(:,:,c)))));
end
</code></pre>

Here’s the result

![image](https://user-images.githubusercontent.com/45420635/49144338-7cad7880-f340-11e8-8e9b-892d23a2c272.png)

![image](https://user-images.githubusercontent.com/45420635/49144347-846d1d00-f340-11e8-914d-db23af5469a0.png)

![image](https://user-images.githubusercontent.com/45420635/49144359-8a62fe00-f340-11e8-8a66-f1314917dd5c.png)

![image](https://user-images.githubusercontent.com/45420635/49144372-8fc04880-f340-11e8-8d53-de273524ef71.png)

In the four results, the logarithmic merging with the tent weighting function method was most probable. 

### Tonemapping
-------------
We used tonemapping for the above logarithmic merging with the tent weighting function method. We used photographic tonemapping and bilateral filtering tonemapping as the method.
* Photographic tonemapping
The photographic tonemapping equation is as follows. I have tried both processing each channel and changing the domain to process only Luminance. As a result, processing each channel was much better.

![image](https://user-images.githubusercontent.com/45420635/49144423-b1213480-f340-11e8-9479-63d78d58446f.png)

<pre><code> % Photographic Tonemapping

I = final_logm_tent;
K = 0.15;
B = 0.95;
eps = 1e-8;

for c=1:3

    I_white = B*max(max(I(:,:,c)));
    I_m = exp(sum(sum(log(I(:,:,c)+eps)))./(M*N));
    I_hat = K.*I(:,:,c)./I_m;

    final_I(:,:,c) = (I_hat.*(1+(I_hat./(I_white)^2)))./(1+I_hat);
    final_I(:,:,c) = 255.*(final_I(:,:,c)./(max(max(final_I(:,:,c)))));
end
</code></pre>

![image](https://user-images.githubusercontent.com/45420635/49144447-c302d780-f340-11e8-9617-8aea45cc9149.png)

* Bilateral Filtering tonemapping
The Bilateral Filtering tonemapping equation is as follows. I have tried both processing each channel and changing the domain to process only Luminance. As a result, processing each channel was much better.
+ Compute the log intensity
+ Compute the base layer using bilatering filtering
+ Compute the detail layer
+ Apply an offset and a scale S to the base
+ Reconstructed the intensity

<pre><code> S = 4;
I = final_logm_tent;

for c=1:3
    I_L = log(I(:,:,c)+eps);
    I_L = (I_L - min(min(I_L)))./(max(max(I_L))-min(min(I_L)));
    I_B = bfilter2(I_L,3);
    I_D = I_L - I_B;
    B_hat = S.*(I_B - max(max(I_B)));
    I_final(:,:,c) = exp(B_hat+I_D);
    I_final(:,:,c) = 255.*(I_final(:,:,c) - min(min(I_final(:,:,c))))./(max(max(I_final(:,:,c)))-min(min(I_final(:,:,c))));
end
</code></pre>

![image](https://user-images.githubusercontent.com/45420635/49144514-ea59a480-f340-11e8-97d9-69f453d41b91.png)

After tonemapping, I still did not get better results than I thought. I thought about those results, the dataset taken when the shutter speed is small (the space where the Toy story frame is placed) is larger than the dataset taken when the shutter speed is high (space with the shelves). Therefore, it seems that only the bright space with light is well processed.

I apologize for the delay in submission due to the conference attendance. I will use all the remaining free late days. Assignment 4 was the most interesting CP assignment.
