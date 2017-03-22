#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

//Screen information
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

//Camera information
vec3 cameraPos(0, 0, -3.001);
mat3 cameraRot(vec3(1,0,0), vec3(0,1,0), vec3(0,0,1));
float yaw = 0;
float focalLength = 500;
float posDelta = 0.01;
float rotDelta = 0.01;

//Scene information
vector<Triangle> triangles;

//The color of the current triangle being drawn
vec3 currentColor;

//The depth buffer for the screen
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

//Floating point inaccuracy constant
float epsilon = 0.00001;

//Data structure holding information for each pixel
struct Pixel {
	int x;
	int y;
	float zinv;
	vec3 pos3d;
};

struct Vertex {
	vec3 position;
};

vec3 lightPos(0,-0.5, -0.7);
vec3 lightPower = 14.0f * vec3(1,1,1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1,1,1);

vec3 currentNormal;
vec3 currentReflectance;
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader(const vec3& v, ivec2& p);
//void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void updateRotationMatrix();

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels);
void DrawPolygonRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShader(const Vertex& v, Pixel& p);
void DrawPolygon_depth(const vector<Vertex>& vertices);

void PixelShader(const Pixel& p);

int main( int argc, char* argv[] )
{
	//Load in the scene
	LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	/*vector<Pixel> vertexPixels(3);
	vertexPixels[0] = {10, 5, 1, vec3(0,0,0)};
	vertexPixels[1] = {5, 10, 0, vec3(1,0,0)};
	vertexPixels[2] = {15, 15, 3, vec3(1,1,1)};
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for(unsigned int row = 0; row < leftPixels.size(); row++)
	{
		printf("Start: (%d, %d, %lf, (%lf, %lf %lf)). End: (%d, %d, %lf, (%lf, %lf %lf)).\n", leftPixels[row].x, leftPixels[row].y, leftPixels[row].zinv, leftPixels[row].illumination.x, leftPixels[row].illumination.y, leftPixels[row].illumination.z, rightPixels[row].x, rightPixels[row].y, rightPixels[row].zinv, rightPixels[row].illumination.x, rightPixels[row].illumination.y, rightPixels[row].illumination.z);	
	}*/

	/*ivec2 a(499,499);
	ivec2 b(125,374);
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a, b, line);
	for(int i = 0; i < line.size(); i++) {
		printf("(%d, %d)\n", line[i].x, line[i].y);
	}
	return 0;*/
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	//get key presses and update camera position
	Uint8* keystate = SDL_GetKeyState(0);
	if(keystate[SDLK_UP])
	{
		//calculate the z-axis vector and move camera in the positive direction
		vec3 forward(cameraRot[2][0], cameraRot[2][1], cameraRot[2][2]);
		cameraPos += posDelta * forward;
	}
	if(keystate[SDLK_DOWN])
	{	
		//calculate the z-axis vector and move camera in the negative direction
		vec3 forward(cameraRot[2][0], cameraRot[2][1], cameraRot[2][2]);
		cameraPos -= posDelta * forward;
		
	}
	if(keystate[SDLK_LEFT])
	{
		//decrease the rotation angle and update rotation matrix 
		yaw -= rotDelta;
		updateRotationMatrix();
	}
	if(keystate[SDLK_RIGHT])
	{
		//increase the rotation angle and update the rotation matrix
		yaw += rotDelta;
		updateRotationMatrix();
	}
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	//reset the depth buffer
	for(int y = 0; y < SCREEN_HEIGHT; y++) {
		for(int x = 0; x < SCREEN_WIDTH; x++) {
			depthBuffer[y][x] = 0;
		}
	}

	//iterate through all of the triangles
	for(unsigned int i = 0; i < triangles.size(); ++i)
	{
		//Update the color to draw with
		currentColor = triangles[i].color;

		//Create a list of the triangles vertices
		vector<Vertex> vertices(3);
		//initalise the vertexes
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		//draw the triangle
		DrawPolygon_depth(vertices);

	}


	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

//Project a 3D poin v to a 2D point p in the image plane
void VertexShader(const vec3& v, ivec2& p)
{
	//Get the relative position of the point from the camera in the current coordinate system 
	vec3 P = (v - cameraPos) * cameraRot;
	//Calculate the 2D coordinates of the point
	p.x = (int) (focalLength * (P.x / P.z) + ((float) SCREEN_WIDTH / 2.0f));
	p.y = (int) (focalLength * (P.y / P.z) + ((float) SCREEN_HEIGHT / 2.0f));
}

//Interpolate all points on a line between points a and b
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	//Calculate the steps that need to be taken in the X and Y direction
	int N = result.size();
	float stepX = (float) (b.x - a.x) / float(max(N-1,1));
	float stepY = (float) (b.y - a.y) / float(max(N-1,1));
	//Set the accumulator vector to the vector a
	vec2 current(a);
	for(int i = 0; i < N; i++)
	{
		//Store the interpolated point in the output result vector
		result[i].x = round(current.x);
		result[i].y = round(current.y);
		//Update the accumulator vector
		current.x = current.x + stepX;
		current.y = current.y + stepY;
	}
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color)
{
	//Compute the number of pixels to draw
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	
	//Interpolate the points on the line between a and b	
	vector<ivec2> line(pixels);
	Interpolate(a,b,line);

	//Draw each pixel on the line
	for(unsigned int i = 0; i < line.size(); ++i) {
		PutPixelSDL(surface, line[i].x, line[i].y, color);
	}
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();

	//Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices(V);
	for(int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}

	//Loop over all vertices and draw the edge from it to the next vertex:
	for(int i = 0; i < V; ++i)
	{
		//The next vector
		int j = (i+1)%V;
		//Draw white line between current vector i and next vector j
		vec3 color(1,1,1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

//Update the cameras coordinate system
void updateRotationMatrix() {
	//Calcuate new columns for the camera's rotation matrix
	cameraRot[0] = vec3(cos(yaw), 0, -sin(yaw));
	cameraRot[2] = vec3(sin(yaw), 0, cos(yaw));
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels)
{
	//maximum and minimum possible values
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();

	//find the minimum and maximum y values of the vertices in the polygon
	for(unsigned int i = 0; i < vertexPixels.size(); i++) 
	{
		if(vertexPixels[i].y < minY)
			minY = vertexPixels[i].y;
		if(vertexPixels[i].y > maxY)
			maxY = vertexPixels[i].y;
	}

	//resize the leftPixels and rightPixels vector lists to the right number of elements
	unsigned int numRows = maxY - minY + 1;
	leftPixels.resize(numRows);
	rightPixels.resize(numRows);

	//set the leftPixels elecment to the max value and the rightPixels elements to the min value
	for(unsigned int i = 0; i < numRows; i++)
	{
		leftPixels[i].x = numeric_limits<int>::max();
		rightPixels[i].x = numeric_limits<int>::min();
	}

	//Fill the leftPixels and rightPixels lists with the appropriate values
	int V = vertexPixels.size();
	for(int i = 0; i < V; ++i)
	{
		//i is the index of the current vertex and j is the index of the next vertex
		int j = (i+1)%V;
		vec3 color(1,1,1);

		//a is the curent vertex and b is the next vertex
		ivec2 a = vertexPixels[i];
		ivec2 b = vertexPixels[j];

		//calculate the number of pixels needed for the line
		ivec2 delta = glm::abs(a - b);
		int pixels = glm::max(delta.x, delta.y) + 1;
		vector<ivec2> line(pixels);

		//Interpolate the points on the line between vertex a and b
		Interpolate(a,b,line);

		//Loop through all points on the line and update the leftPixels and rightPixels lists
		for(unsigned int k = 0; k < line.size(); k++)
		{
			//if the x coordinate for this y is less than the current x coordinate for this y in 
			//the leftPixels array, then update that x coordinate
			if(leftPixels[line[k].y - minY].x > line[k].x)
			{
				leftPixels[line[k].y - minY].x = line[k].x;
				leftPixels[line[k].y - minY].y = line[k].y;
			}
			//if the x coordinate for this y is greater than the current x coordinate for this y in
			//the rightPixels array, then update that x coordinate
			if(rightPixels[line[k].y - minY].x < line[k].x)
			{
				rightPixels[line[k].y - minY].x = line[k].x;
				rightPixels[line[k].y - minY].y = line[k].y;
			}
		}
	}
}

//given the left most and right most x coordinates for the corresponding y coordinate, draw the polygon
void DrawPolygonRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels)
{
	for(unsigned int i = 0; i < leftPixels.size(); i++)
	{
		for(int x = leftPixels[i].x; x <= rightPixels[i].x; x++)
		{
			//Check that we arent drawing outside of the screen
			if(x >= 0 && x <= SCREEN_WIDTH && leftPixels[i].y >= 0 && leftPixels[i].y <= SCREEN_HEIGHT)
				PutPixelSDL(screen, x, leftPixels[i].y, currentColor);
		}
	}
}


//Draw a polygon given its vertices
void DrawPolygon(const vector<vec3>& vertices)
{
	//create a list to store the projected 2D positions of the vertexes
	int V = vertices.size();
	vector<ivec2> vertexPixels(V);

	//for each vertex calculate its projected 2D coordinate
	for(int i = 0; i < V; i++) {
		VertexShader(vertices[i], vertexPixels[i]);
	}	

	//calculate leftmost and rightmost pixels of the polygon for each y position
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);

	//draw the polygon using these pixels coordinates
	DrawPolygonRows(leftPixels, rightPixels);

}


















//Interpolate the points on a line between vectors a and b and put the points into result
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {

	//Calculate the steps needed in the x and y direction, and also the step needed for the zinv
	int N = result.size();
	float stepX = (float) (b.x - a.x) / float(max(N-1,1));
	float stepY = (float) (b.y - a.y) / float(max(N-1,1));
	float stepZinv = (b.zinv - a.zinv) / float(max(N-1,1));
	vec3 stepPos3d = (b.pos3d - a.pos3d) / float(max(N-1,1));

	//Set the accumulator values to the same as vector a
	float currentX = (float) a.x;
	float currentY = (float) a.y;
	float currentZinv = (float) a.zinv;
	vec3 currentPos3d = a.pos3d;

	//Calculate each point on the line
	for(int i = 0; i < N; i++)
	{
		//Store the interpolated point in the result list
		result[i].x = round(currentX);
		result[i].y = round(currentY);
		result[i].zinv = currentZinv;
		result[i].pos3d = currentPos3d;

		//Update the accumulator values
		currentX += stepX;
		currentY += stepY;
		currentZinv += stepZinv;
		currentPos3d += stepPos3d;
	}
}

//Calculate the leftPixels and rightPixels lists from the vertex coordinats of the polygon
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
		
	//The maximum and minimum integer values
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();

	//Calculate the smallest and largest y values of the polygon
	for(unsigned int i = 0; i < vertexPixels.size(); i++) 
	{
		if(vertexPixels[i].y < minY)
			minY = vertexPixels[i].y;
		if(vertexPixels[i].y > maxY)
			maxY = vertexPixels[i].y;
	}

	//resize the leftPixels and rightPixels lists to the correct size
	unsigned int numRows = maxY - minY + 1;
	leftPixels.resize(numRows);
	rightPixels.resize(numRows);

	//set the leftPixels x values to the maximum integer value and
	//the rightPixels x values to the minimum integer value
	for(unsigned int i = 0; i < numRows; i++)
	{
		leftPixels[i].x = numeric_limits<int>::max();
		rightPixels[i].x = numeric_limits<int>::min();
	}


	//Iterate through all of the vertexes
	int V = vertexPixels.size();
	for(int i = 0; i < V; ++i)
	{
		//i is the index of the current vertex and j is the index of the next vertex
		int j = (i+1)%V;
		
		//a is the current vertex and b is the next vertex
		Pixel a = vertexPixels[i];
		Pixel b = vertexPixels[j];

		//calculate the number of pixels needed for the line between a and b		
		int pixels = max(abs(a.x - b.x), abs(a.y - b.y)) + 1;
		vector<Pixel> line(pixels);

		//Interpolate the points on the line between a and b
		Interpolate(a,b,line);
	
		//Iterate through all of the interpolated points on the line		
		for(unsigned int k = 0; k < line.size(); k++)
		{

			//if the x value for this y coordinate is less than the current x value for
			//this y coordinate then update the entry with the current point values
			if(leftPixels[line[k].y - minY].x > line[k].x)
			{
				leftPixels[line[k].y - minY].x = line[k].x;
				leftPixels[line[k].y - minY].y = line[k].y;
				leftPixels[line[k].y - minY].zinv = line[k].zinv;
				leftPixels[line[k].y - minY].pos3d = line[k].pos3d;
			}

			//if the x value for this y coordinate is greater than the current x value for
			//this y coordinate then update the entry with the current point values
			if(rightPixels[line[k].y - minY].x < line[k].x)
			{
				rightPixels[line[k].y - minY].x = line[k].x;
				rightPixels[line[k].y - minY].y = line[k].y;
				rightPixels[line[k].y - minY].zinv = line[k].zinv;
				rightPixels[line[k].y - minY].pos3d = line[k].pos3d;
			}
		}
	}
}

//Draw the polygon given the leftmost and rightmost pixel positions for each y value
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
	//Iterate through all the leftPixels and rightPixel elements
	unsigned int V = leftPixels.size();
	for(unsigned int i = 0; i < V; i++) {

		//Calculate the number of pixels needed for the line between corresponding leftPixels and rightPixels elements
		int pixels = rightPixels[i].x - leftPixels[i].x + 1;
		vector<Pixel> line(pixels);

		//Interpolate between these pixels (as we need to calculate the zinv values for each pixel)
		Interpolate(leftPixels[i], rightPixels[i], line);

		//Iterate over all the interpolated pixels
		for(int j = 0; j < pixels; j++) {
			//If the pixels zinv value is greater than the corresponding one in the buffer, then draw the pixel
			PixelShader(line[j]);
		}
	}
}

//Calculate the Euclidean distance between the two given vectors
float distanceBetweenPoints(vec3 a, vec3 b) {
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

//Caluclate the dot product between the two given vectors
float dotProduct(vec3 a, vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

//Calculate the normalised direction vector from the given vector to the lightsource
vec3 unitVectorToLightSource(vec3 a) {
	vec3 v(lightPos.x - a.x, lightPos.y - a.y, lightPos.z - a.z);
	return normalize(v);
}


//Project the 3d point v to 2D
void VertexShader(const Vertex& v, Pixel& p) {

	//get the relative position of v from C in the current coordinate system
	vec3 C = (v.position - cameraPos) * cameraRot;

	//project the x and y values
	p.x = (int) (focalLength * (C.x / C.z) + ((float) SCREEN_WIDTH / 2.0f));
	p.y = (int) (focalLength * (C.y / C.z) + ((float) SCREEN_HEIGHT / 2.0f));

	//calculate the inverse of the depth of the point
	p.zinv = 1.0f/(float)C.z;

	p.pos3d = v.position * p.zinv;
}

void PixelShader(const Pixel& p) {
	
	//Pi constant
	const float pi = 3.1415926535897;

	vec3 pos3d = p.pos3d / p.zinv;

	//distance from intersection point to light source
	float radius = distanceBetweenPoints(pos3d, lightPos);

	//The power per area at this point
	vec3 B = lightPower / (4 * pi * pow(radius,3));

	//unit vector describing normal of surface
	vec3 n = currentNormal;

	//unit vector describing direction from surface point to light source
	vec3 r = unitVectorToLightSource(pos3d);

	//fraction of the power per area depending on surface's angle from light source
	vec3 D = B * max(dotProduct(r,n),0.0f);

	vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);

	//If pixels depth is less than the current pixels depth in the image
	//then update the image
	if(p.zinv > depthBuffer[p.y][p.x] + epsilon) {
		depthBuffer[p.y][p.x] = p.zinv;
		PutPixelSDL(screen, p.x, p.y, illumination);
	}
}

//Draw a polygon given its vertices, taking into account depth
void DrawPolygon_depth(const vector<Vertex>& vertices)
{

	//Calculate the projection of the polygons vertexes
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for(int i = 0; i < V; i++) {
		VertexShader(vertices[i], vertexPixels[i]);
	}
	
	//lists to store the left most and right post pixel x values for each y value
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;

	//Compute the leftPixels and rightPixels lists from the vertexes
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);

	//Draw the polygon
	DrawPolygonRows(leftPixels, rightPixels);
}