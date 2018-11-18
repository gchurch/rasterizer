#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

class V
{
public:
	glm::vec3 pos3d;
	glm::ivec2 textureCoordinates;

	V(glm::vec3 pos3d, glm::ivec2 textureCoordinates)
		: pos3d(pos3d), textureCoordinates(textureCoordinates)
	{

	}
};

// Used to describe a triangular surface:
class Triangle
{
public:
	V v0;
	V v1;
	V v2;
	glm::vec3 normal;
	glm::vec3 color;
	int texture;

	Triangle( V v0, V v1, V v2, glm::vec3 color, int texture )
		: v0(v0), v1(v1), v2(v2), color(color), texture(texture)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
		glm::vec3 e1 = v1.pos3d-v0.pos3d;
		glm::vec3 e2 = v2.pos3d-v0.pos3d;
		normal = glm::normalize( glm::cross( e2, e1 ) );
	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
	using glm::vec3;
	using glm::ivec2;

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );

	triangles.clear();
	triangles.reserve( 5*2*3 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	V A(vec3(L,0,0), ivec2(511,511));
	V B(vec3(0,0,0), ivec2(0,511));
	V C(vec3(L,0,L), ivec2(511,0));
	V D(vec3(0,0,L), ivec2(0,0));

	V E(vec3(L,L,0), ivec2(0,0));
	V F(vec3(0,L,0), ivec2(0,0));
	V G(vec3(L,L,L), ivec2(0,0));
	V H(vec3(0,L,L), ivec2(0,0));


	// Floor:
	triangles.push_back( Triangle( C, B, A, green, 1 ) );
	triangles.push_back( Triangle( C, D, B, green, 1 ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, purple, 0 ) );
	triangles.push_back( Triangle( C, E, G, purple, 0 ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, yellow, 0 ) );
	triangles.push_back( Triangle( H, F, D, yellow, 0 ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, cyan, 0 ) );
	triangles.push_back( Triangle( F, H, G, cyan, 0 ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, white, 0 ) );
	triangles.push_back( Triangle( G, H, D, white, 0 ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = V(vec3(290,0,114), ivec2(0,0));
	B = V(vec3(130,0, 65), ivec2(0,0));
	C = V(vec3(240,0,272), ivec2(0,0));
	D = V(vec3( 82,0,225), ivec2(0,0));

	E = V(vec3(290,165,114), ivec2(0,0));
	F = V(vec3(130,165, 65), ivec2(0,0));
	G = V(vec3(240,165,272), ivec2(0,0));
	H = V(vec3( 82,165,225), ivec2(0,0));

	// Front
	triangles.push_back( Triangle(E,B,A,red,0) );
	triangles.push_back( Triangle(E,F,B,red,0) );

	// Front
	triangles.push_back( Triangle(F,D,B,red,0) );
	triangles.push_back( Triangle(F,H,D,red,0) );

	// BACK
	triangles.push_back( Triangle(H,C,D,red,0) );
	triangles.push_back( Triangle(H,G,C,red,0) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,red,0) );
	triangles.push_back( Triangle(E,A,C,red,0) );

	// TOP
	triangles.push_back( Triangle(G,F,E,red,0) );
	triangles.push_back( Triangle(G,H,F,red,0) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = V(vec3(423,0,247), ivec2(0,0));
	B = V(vec3(265,0,296), ivec2(0,0));
	C = V(vec3(472,0,406), ivec2(0,0));
	D = V(vec3(314,0,456), ivec2(0,0));

	E = V(vec3(423,330,247), ivec2(0,0));
	F = V(vec3(265,330,296), ivec2(0,0));
	G = V(vec3(472,330,406), ivec2(0,0));
	H = V(vec3(314,330,456), ivec2(0,0));

	// Front
	triangles.push_back( Triangle(E,B,A,blue,0) );
	triangles.push_back( Triangle(E,F,B,blue,0) );

	// Front
	triangles.push_back( Triangle(F,D,B,blue,0) );
	triangles.push_back( Triangle(F,H,D,blue,0) );

	// BACK
	triangles.push_back( Triangle(H,C,D,blue,0) );
	triangles.push_back( Triangle(H,G,C,blue,0) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,blue,0) );
	triangles.push_back( Triangle(E,A,C,blue,0) );

	// TOP
	triangles.push_back( Triangle(G,F,E,blue,0) );
	triangles.push_back( Triangle(G,H,F,blue,0) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0.pos3d *= 2/L;
		triangles[i].v1.pos3d *= 2/L;
		triangles[i].v2.pos3d *= 2/L;

		triangles[i].v0.pos3d -= vec3(1,1,1);
		triangles[i].v1.pos3d -= vec3(1,1,1);
		triangles[i].v2.pos3d -= vec3(1,1,1);

		triangles[i].v0.pos3d.x *= -1;
		triangles[i].v1.pos3d.x *= -1;
		triangles[i].v2.pos3d.x *= -1;

		triangles[i].v0.pos3d.y *= -1;
		triangles[i].v1.pos3d.y *= -1;
		triangles[i].v2.pos3d.y *= -1;

		triangles[i].ComputeNormal();
	}
}

#endif