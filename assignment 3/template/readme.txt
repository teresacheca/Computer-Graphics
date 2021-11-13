Authors: Teresa del Carmen Checa Marabotto
	 Lucía Sánchez-Montes Gómez

EXERCISE 1

In order to do this exercise we need to use the Ray-plane intersection formula, we didn’t have any problem with this part. But it was important to remember that we needed to initialize the values.
In the creation of the planes we have more problems because at first we didn’t know how to find the values, when we found them we had the problem that when we printed the planes they were double. The solution to this was checking the value of t and if it was minus from 0 there was no hit.


EXERCISE 2

checkerboardTexture and rainbowTexture


To solve this exercise, we first had to calculate the uv coordinates of the texture. We found some problems doing this, mainly between the parameters for asin and atan2 in terms of degrees and radians, but figured out that our result did not have to be converted, since it was as we wanted. Then, we had to apply the clamp function to make sure that the variables were within the range of values they should.


After this, we had to modify the values of the material texture, to modify the phongModel. The biggest problem with this part was understanding how we should call the texture of a material and also that it should be added to the color of the figure.


The last part of the exercise consisted of modifying the file texture.h
We solved both the checkerboardTexture and the rainbowTexture.
In this part, we had to apply the formulas correctly and vary their parameters to obtain the expected results (more or less squares or lines and their inclinations) as well as the modification of the colors in the rainbowTExture to get the colors we wanted.




EXERCISE 3
There are two parts in this exercise. The first one is to extend the Phong Model formula by adding the attenuation of the light. The second part is the tonemapped in this part we have to apply the formula given in the pdf. In this exercise we have a with the pow function because at first I was trying to power the whole vector but I couldn’t so I did it coordinate by coordinate and it worked.
