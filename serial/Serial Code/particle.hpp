#include "common.hpp"

class particle
{
public:
	float mass;
	double velocityX;
	double velocityY;
	double xPosn;
	double yPosn;
	double force;
	double forceX;
	double forceY;
	bool hasLeftSystem;
	particle();
	void generateXPosn();
	void generateYPosn();
	void updateParticle();
};


