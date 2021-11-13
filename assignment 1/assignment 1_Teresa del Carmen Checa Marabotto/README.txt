AUTHOR: TERESA DEL CARMEN CHECA MARABOTTO

EXERCISES SOLVED: 1, 2 AND 3

EXERCISE 1:
	Author: Teresa del Carmen Checa Marabotto

	Encountered problems: non. Just need to remember the correct formula



EXERCISE 2:
	Author: Teresa del Carmen Checa Marabotto

	Encountered problems: 
	
		- Do the normal of c: at first I was using the function called "glm::normalize()" but I realize that it returned a vector, not a number, so I decided to put the exact formula for that: sqrt( pow(c[0],2) + pow(c[1],2) + pow(c[2],2) );

		- Calculate the distance: at first I was using the normal of gamma function to calculate the distance between the origin of the ray to the intersection point, but I discovered that I didn't have snese because if we use the normal, we are not calculating a distance, so I used just the gamma function to calculate it.

		- Working with any origin: at first I was using the origin (0,0,0) to do the calculations, but later I realized that it need to be changed for ray.origin, so it can work for others types of origins diferents from (0,0,0)

		- Start hit: it's take me so long to remember that hit.hit had to be changed to true

 

EXERCISE 3:

	Author: Teresa del Carmen Checa Marabotto

	Encountered problems:
		
		- Intersection part: The intersection part was not calculated correctly, because I was not comparing distances. I had to calculate the distances of each point separately. Then, check which one is the closest to the origin, that is, the smallest distance, and modify hit according to that point

