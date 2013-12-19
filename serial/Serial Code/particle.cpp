#include "particle.hpp"

particle::particle()
{
	mass = MASS;
	velocityX = INITVELOCITY;
	velocityY = INITVELOCITY;
	float xPosn = XMIN;
	float yPosn = YMIN;
	force = 0.0;
	forceX = 0.0;
	forceY = 0.0;
	hasLeftSystem = false;
}

void particle::generateXPosn()
{
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = XDIM - XMIN;
	float r = random * diff;
	xPosn = r;
}

void particle::generateYPosn()
{
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = YDIM - YMIN;
	float r = random * diff;
	yPosn = r;
}

void particle::updateParticle()
{
	//update velocities Using v = u + at
	velocityX = velocityX + ((forceX / MASS)*TIMESTEP);
	velocityY = velocityY + ((forceY / MASS)*TIMESTEP);
	//update positions
	xPosn += (velocityX*TIMESTEP);
	yPosn += (velocityY*TIMESTEP);
		
	if (xPosn >= XDIM || xPosn < XMIN)
	{
		//cout << "\n particle went out of the x dimensions :(";
		hasLeftSystem = true;
	}
	if (yPosn >= YDIM || yPosn < YMIN)
	{
		//cout << "\n particle went out of the Y dimensions :(";
		hasLeftSystem = true;
	}
}

