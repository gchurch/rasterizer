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

vec3 currentColor;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader(const vec3& v, ivec2& p);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void updateRotationMatrix();
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels);
void DrawPolygonRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);

int main( int argc, char* argv[] )
{
	//Load in the scene
	/*LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );*/

	vector<ivec2> vertexPixels(3);
	vertexPixels[0] = ivec2(374, 374);
	vertexPixels[1] = ivec2(374, 374);
	vertexPixels[2] = ivec2(500, 500);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for(unsigned int row = 0; row < leftPixels.size(); row++)
	{
		cout << "Start: ("
		<< leftPixels[row].x << ","
		<< leftPixels[row].y << "). "
		<< "End: ("
		<< rightPixels[row].x << ","
		<< rightPixels[row].y << "). " << endl;
	}
	return 0;
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

	for(unsigned int i = 0; i < triangles.size(); ++i)
	{
		currentColor = triangles[i].color;

		vector<vec3> vertices(3);
		
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		DrawPolygon(vertices);

	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const vec3& v, ivec2& p)
{
	vec3 P = (v - cameraPos) * cameraRot;
	p.x = (int) (focalLength * (P.x / P.z) + ((float) SCREEN_WIDTH / 2.0f));
	p.y = (int) (focalLength * (P.y / P.z) + ((float) SCREEN_HEIGHT / 2.0f));
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current(a);
	for(int i = 0; i < N; ++i)
	{
		result[i] = current;
		current += step;
	}
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color)
{
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a,b,line);
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
		int j = (i+1)%V; //the next vertex
		vec3 color(1,1,1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

void updateRotationMatrix() {
	//Calcuate new columns for the camera's rotation matrix
	cameraRot[0] = vec3(cos(yaw), 0, -sin(yaw));
	cameraRot[2] = vec3(sin(yaw), 0, cos(yaw));
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels)
{
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();

	for(unsigned int i = 0; i < vertexPixels.size(); i++) 
	{
		printf("v%d: (%d, %d), ", i, vertexPixels[i].x, vertexPixels[i].y);

		if(vertexPixels[i].y < minY)
			minY = vertexPixels[i].y;
		if(vertexPixels[i].y > maxY)
			maxY = vertexPixels[i].y;
	}
	printf("\n");

	leftPixels.resize(maxY - minY + 1);
	rightPixels.resize(maxY - minY + 1);

	for(unsigned int i = 0; i < leftPixels.size(); i++)
	{
		leftPixels[i].x = numeric_limits<int>::max();
		rightPixels[i].x = numeric_limits<int>::min();
	}

	int V = vertexPixels.size();

	for(int i = 0; i < V; ++i)
	{
		int j = (i+1)%V;
		vec3 color(1,1,1);
		ivec2 a = vertexPixels[i];
		ivec2 b = vertexPixels[j];
		ivec2 delta = glm::abs(a - b);
		int pixels = glm::max(delta.x, delta.y) + 1;
		vector<ivec2> line(pixels);
		Interpolate(a,b,line);
		for(unsigned int k = 0; k < line.size(); k++)
		{
			printf("%d ", k);
			if(leftPixels[line[k].y - minY].x > line[k].x) 
			{
				leftPixels[line[k].y - minY].x = line[k].x; 
				leftPixels[line[k].y - minY].y = line[k].y;
			}
			if(rightPixels[line[k].y - minY].x < line[k].x)
			{
				rightPixels[line[k].y - minY].x = line[k].x;
				rightPixels[line[k].y - minY].y = line[k].y;
			}
			printf("done\n");
		}
		printf("fin\n");
	}

	printf("finished\n");
}

void DrawPolygonRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels)
{
	for(unsigned int i = 0; i < leftPixels.size(); i++)
	{
		for(int x = leftPixels[i].x; x <= rightPixels[i].x; x++)
		{
			if(x >= 0 && x <= SCREEN_WIDTH && leftPixels[i].y >= 0 && leftPixels[i].y <= SCREEN_HEIGHT)
				PutPixelSDL(screen, x, leftPixels[i].y, currentColor);
		}
	}
}

void DrawPolygon(const vector<vec3>& vertices)
{
	int V = vertices.size();

	vector<ivec2> vertexPixels(V);
	for(int i = 0; i < V; i++)
		VertexShader(vertices[i], vertexPixels[i]);
	
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);

}
