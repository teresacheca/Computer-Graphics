/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

#include "Image.h"
#include "Material.h"

using namespace std;

/**
 Class representing a single ray.
 */
class Ray{
public:
    glm::vec3 origin; ///< Origin of the ray
    glm::vec3 direction; ///< Direction of the ray
	/**
	 Contructor of the ray
	 @param origin Origin of the ray
	 @param direction Direction of the ray
	 */
    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction){
    }
};


class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit{
    bool hit; ///< Boolean indicating whether there was or there was no intersection with an object
    glm::vec3 normal; ///< Normal vector of the intersected object at the intersection point
    glm::vec3 intersection; ///< Point of Intersection
    float distance; ///< Distance from the origin of the ray to the intersection point
    Object *object; ///< A pointer to the intersected object
	glm::vec2 uv; ///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object{
	
protected:
	glm::mat4 transformationMatrix; ///< Matrix representing the transformation from the local to the global coordinate system
	glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
	glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the local to the global coordinate system
	
public:
	glm::vec3 color; ///< Color of the object
	Material material; ///< Structure describing the material of the object
	/** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

	/** Function that returns the material struct of the object*/
	Material getMaterial(){
		return material;
	}
	/** Function that set the material
	 @param material A structure describing the material of the object
	*/
	void setMaterial(Material material){
		this->material = material;
	}
	/** Functions for setting up all the transformation matrices
	 @param matrix The matrix representing the transformation of the object in the global coordinates */
	void setTransformation(glm::mat4 matrix){
		
		transformationMatrix = matrix;
		
		
		/* ----- Assignment 5 ---------
		 Set the two remaining matrices
		
		inverseTransformationMatrix =
		normalMatrix =
		 
		 */
		inverseTransformationMatrix = glm::inverse(matrix);
		normalMatrix = glm::transpose(inverseTransformationMatrix);
	}
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object{
private:
    float radius; ///< Radius of the sphere
    glm::vec3 center; ///< Center of the sphere

public:
	/**
	 The constructor of the sphere
	 @param radius Radius of the sphere
	 @param center Center of the sphere
	 @param color Color of the sphere
	 */
    Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center){
		this->color = color;
    }
	Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center){
		this->material = material;
	}
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){

        glm::vec3 c = center - ray.origin;

        float cdotc = glm::dot(c,c);
        float cdotd = glm::dot(c, ray.direction);

        Hit hit;

        float D = 0;
		if (cdotc > cdotd*cdotd){
			D =  sqrt(cdotc - cdotd*cdotd);
		}
        if(D<=radius){
            hit.hit = true;
            float t1 = cdotd - sqrt(radius*radius - D*D);
            float t2 = cdotd + sqrt(radius*radius - D*D);

            float t = t1;
            if(t<0) t = t2;
            if(t<0){
                hit.hit = false;
                return hit;
            }

			hit.intersection = ray.origin + t * ray.direction;
			hit.normal = glm::normalize(hit.intersection - center);
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
			
			hit.uv.s = (asin(hit.normal.y) + M_PI/2)/M_PI;
			hit.uv.t = (atan2(hit.normal.z,hit.normal.x) + M_PI) / (2*M_PI);
        }
		else{
            hit.hit = false;
		}
		return hit;
    }
};


class Plane : public Object{

private:
	glm::vec3 normal;
	glm::vec3 point;

public:
	Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal){
	}
	Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;
		float DdotN = glm::dot(ray.direction, normal);
		if(DdotN < 0){
			
			float PdotN = glm::dot (point-ray.origin, normal);
			float t = PdotN/DdotN;
			
			if(t > 0){
				hit.hit = true;
				hit.normal = normal;
				hit.distance = t;
				hit.object = this;
				hit.intersection = t * ray.direction + ray.origin;
			}
		}
		return hit;
	}
};

class Cone : public Object{
public:
	Cone(Material material){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;
		
	
		/*  ---- Assignemnt 5 -----
		
		 Implement the ray-cone intersection. Before intersecting the ray with the cone,
		 make sure that you transform the ray into the local coordinate system.
		 Remember about normalizing all the directions after transformations.
		 
		*/
	
		/* If the intersection is found, you have to set all the critical fields in the Hit strucutre
		 Remember that the final information about intersection point, normal vector and distance have to be given
		 in the global coordinate system.
		 */

		hit.hit = true;
		hit.object = this;

		//The plane that will use to create the disk that closes the cone:
		Plane close = Plane(glm::vec3(0.f,1.f,0.f),glm::vec3(0.f, 1.f,0.f),material);

		ray.direction = glm::normalize(ray.direction);

		//transform to local coordinates
		glm::vec3 newOrigin = inverseTransformationMatrix * glm::vec4(ray.origin, 1.0);
		glm::vec3 newDirection = inverseTransformationMatrix * glm::vec4(ray.direction, 0.0);
		newDirection = glm::normalize(newDirection);
	
		float a, b, c;
		a = pow(newDirection.x, 2) + pow(newDirection.z, 2) - pow(newDirection.y, 2);
		b = 2 * ((newOrigin.x * newDirection.x) + (newOrigin.z * newDirection.z) - (newOrigin.y * newDirection.y));
		c = pow(newOrigin.x, 2) + pow(newOrigin.z, 2) - pow(newOrigin.y, 2);
	  

		float delta = pow(b, 2) - 4 * a * c;
		float t=0;
		if(delta < 0){
			hit.hit = false;
			return hit;
		}
		if(delta==0){
			t=-b/(2*a);
		}
		else
		{
			float t1 = (-b + sqrt(delta))/ (2*a);
			float t2 = (-b - sqrt(delta)) / (2*a);
			
			t = t1;
            //We will check t1 and t2 and get the smallest, because the other one will be behind and we will not see it
			if(t2<=t1) t=t2;

			//The insertection can not be < 0, because that means that it is 
			// behind us, so, it will not hit the figure
            if(t<0){
                hit.hit = false;
                return hit;
            }

			//Calculation of the new intersection:
			hit.intersection = newOrigin + t * newDirection;
		
			//We will check the the y component is not smaller than 0. If it is, we will not intersect the cone
			//because it will create an infinite cone
			if(hit.intersection.y < 0 ){
				hit.hit = false;
				return hit;
			}

			//EXERCISE 2: creation of the disk which closes the cone
			if( hit.intersection.y > 1){
				//if we are in the edge (we are in local coordinates so is 1) we know that the plane will intersect 
				//with the cone so we make a call to the ray intersect with the plane
				hit = close.intersect(Ray(newOrigin,newDirection));
				if(hit.hit){
					float check = (float) sqrt(pow(hit.intersection.x,2)+pow(hit.intersection.z,2));
					if(check>1){
						hit.hit=false;
						return hit;
					}

					//Change to global coordinates before return it
					hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1.0);
					hit.normal = normalMatrix * glm::vec4(hit.normal,0.0);
					hit.distance = glm::distance(ray.origin, hit.intersection);
					hit.object = this;
					return hit;
				}
			}


			//Calculation of the rest of the parameters of the hit
			hit.normal=glm::vec3(hit.intersection.x,-hit.intersection.y, hit.intersection.z );
			hit.normal = glm::normalize(hit.normal);
			hit.distance = glm::distance(newOrigin, hit.intersection);
		}

		//transform to global coordinates
		hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1.0);
		hit.normal = normalMatrix * glm::vec4(hit.normal,0.0);
		hit.normal = glm::normalize(hit.normal);
		hit.distance = glm::distance(ray.origin, hit.intersection);

		

		return hit;
	}
};


/**
 Light class
 */
class Light{
public:
	glm::vec3 position; ///< Position of the light source
	glm::vec3 color; ///< Color/intentisty of the light source
	Light(glm::vec3 position): position(position){
		color = glm::vec3(1.0);
	}
	Light(glm::vec3 position, glm::vec3 color): position(position), color(color){
	}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001,0.001,0.001);
vector<Object *> objects; ///< A list of all objects in the scene


/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction, Material material){

	glm::vec3 color(0.0);
	for(int light_num = 0; light_num < lights.size(); light_num++){

		glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
		glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

		float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
		float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

		
		glm::vec3 diffuse_color = material.diffuse;
		if(material.texture){
			diffuse_color = material.texture(uv);
		}
		
		glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
		glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));
		
		
		// distance to the light
		float r = glm::distance(point,lights[light_num]->position);
		r = max(r, 0.1f);
		

		color += lights[light_num]->color * (diffuse + specular) / r/r;
	}
	color += ambient_light * material.ambient;
	
	color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
	return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray){

	Hit closest_hit;

	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	for(int k = 0; k<objects.size(); k++){
		Hit hit = objects[k]->intersect(ray);
		if(hit.hit == true && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}

	glm::vec3 color(0.0);

	if(closest_hit.hit){
		color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
	}else{
		color = glm::vec3(0.0, 0.0, 0.0);
	}
	return color;
}
/**
 Function defining the scene
 */
void sceneDefinition (){

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
	green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
	red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;


	objects.push_back(new Sphere(1.0, glm::vec3(1,-2,8), blue_specular));
	objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));
	
	
	
	// ------ Assignment 5 -------
	
	// You can remove the green sphere as it should be replaced with a green cone
	//objects.push_back(new Sphere(1.0, glm::vec3(3,-2,6), green_diffuse));
	
	
	
	//Textured sphere
	Material textured;
	textured.texture = &rainbowTexture;
	objects.push_back(new Sphere(7.0, glm::vec3(-6,4,23), textured));
	
	
	//Planes
	Material red_diffuse;
	red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
	red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);
		
	Material blue_diffuse;
	blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
	blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
	objects.push_back(new Plane(glm::vec3(0,-3,0), glm::vec3(0.0,1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,30), glm::vec3(0.0,0.0,-1.0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(-15,1,0), glm::vec3(1.0,0.0,0.0), red_diffuse));
	objects.push_back(new Plane(glm::vec3(15,1,0), glm::vec3(-1.0,0.0,0.0), blue_diffuse));
	objects.push_back(new Plane(glm::vec3(0,27,0), glm::vec3(0.0,-1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,-0.01), glm::vec3(0.0,0.0,1.0), green_diffuse));

	
	/* ----- Assignment 5 -------
	Create two conse and add them to the collection of our objects.
	Remember to create them with corresponding materials and transformation matrices
	
	
	*/

	//YELLOW CONE:

	//Material for the yellow cone
	Material yellow;
	yellow.ambient = glm::vec3(1.0f, 1.0f, 0.0f);
	yellow.diffuse = glm::vec3(1.0f, 1.0f, 0.0f);
	yellow.specular = glm::vec3(0.7);
	yellow.shininess = 10.0;

	//Transformation matric for the yellow cone
	glm::mat4 m = glm::mat4(1.0f);
	m = glm::translate(m, glm::vec3(5.0f, 9.0f,14.0f));
	m = glm::scale(m, glm::vec3(3.0f, 12.0f, 3.0f));
	m = glm::rotate(m, glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	
	//Creation of the yellow cone
	Cone *cone1 = new Cone(yellow);
	cone1->setTransformation(m);
	objects.push_back(cone1);

	//GREEN CONE:

	//Material for the green cone
	Material green;
	green.ambient = glm::vec3(0.6f, 0.6f, 0.6f);
	green.diffuse = glm::vec3(0.0f, 1.0f, 0.0f);
	green.specular = glm::vec3(0.7);
	green.shininess = 100.0;

	//Transformation matrix for the green cone
	glm::mat4 m2 = glm::mat4(1.0f);
	m2 = glm::translate(m2, glm::vec3(6.0f, -3.0f, 7.0f));
	m2 = glm::rotate(m2, glm::radians(70.0f), glm::vec3(0.0f, 0.0f, 1.0f));
	m2 = glm::scale(m2, glm::vec3(1.0f, 3.0f, 1.0f));

	//Creation of the green cone
	Cone *cone2 = new Cone(green);
	cone2->setTransformation(m2);
	objects.push_back(cone2);



	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity){
	float gamma = 1.0/2.0;
	float alpha = 12.0f;
	return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char * argv[]) {

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; //width of the image
    int height = 768; // height of the image
    float fov = 90; // field of view

	sceneDefinition(); // Let's define a scene

	Image image(width,height); // Create an image where we will store the result

    float s = 2*tan(0.5*fov/180*M_PI)/width;
    float X = -s * width / 2;
    float Y = s * height / 2;

    for(int i = 0; i < width ; i++)
        for(int j = 0; j < height ; j++){

			float dx = X + i*s + s/2;
            float dy = Y - j*s - s/2;
            float dz = 1;

			glm::vec3 origin(0, 0, 0);
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);

            Ray ray(origin, direction);

			image.setPixel(i, j, toneMapping(trace_ray(ray)));

        }

    t = clock() - t;
    cout<<"It took " << ((float)t)/CLOCKS_PER_SEC<< " seconds to render the image."<< endl;
    cout<<"I could render at "<< (float)CLOCKS_PER_SEC/((float)t) << " frames per second."<<endl;

	// Writing the final results of the rendering
	if (argc == 2){
		image.writeImage(argv[2]);
	}else{
		image.writeImage("./result.ppm");
	}

    return 0;
}