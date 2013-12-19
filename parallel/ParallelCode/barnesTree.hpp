#include "node.cpp"
#include <limits>

typedef std::numeric_limits< double > dbl;

class barnesTree
{
public:
	node *m_root;
	long m_noOfParticles;
	long m_noOfSteps;
	float m_theta;
	particle **particleList;
	ofstream logFile;

	barnesTree();
	barnesTree(long noOfParticles, long noOfSteps, float theta);
	void prepareInitList();
	bool run();
	void rebuildTree();
	void clearTree(node *currNode);
	void insert(node *currNode,const int idx);
	void parallelCalcCenterOfMass();
	void calcCenterOfMass(node *currNode);
	void calcForce(node *currNode,long int idx);
	bool checkSDRatio(node *currNode, particle *currParticle);
	int checkQuadrant(node *currNode, particle *temp);
	void splitNode(node *currNode, int idx);
	void logParticlePosition(int stepNumber);
};
