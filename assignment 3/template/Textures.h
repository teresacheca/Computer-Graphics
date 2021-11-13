//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h


#include "glm/glm.hpp"


/*
	Authors: Lucía Sánchez-Montes Gómez
		Teresa del Carmen Checa Marabotto
*/

//We had solved both textures

glm::vec3 checkerboardTexture(glm::vec2 uv){

    /*
     
     
        Exercise 2 (3 points)
     
     
    */
	
	float f =( (int)( floor(30*uv.s) + floor(30*uv.t) ) %2);
	

    return glm::vec3(f);
}
glm::vec3 rainbowTexture(glm::vec2 uv){
    /*
     
     
        Exercise 2 (5 points)
     
     
    */
    glm::vec3 sol;
	
	float f =( (int)( floor(3*uv.s * 30*uv.t) ) %3);

	if( f == 0){
		sol[0] = 1.0;
		sol[1] = 0.0;
		sol[2] = 0.0;

	}else if(f == 1 ){
		sol[0] = 0.0;
		sol[1] = 1.0;
		sol[2] = 0.0;
	}else{
		sol[0] = 0.0;
		sol[1] = 0.0;
		sol[2] = 1.0;
	}

	

    return sol;
}

#endif /* Textures_h */
