#include "particle.cpp"

class node
{
public:
        node *nwChild;
        node *neChild;
        node *swChild;
        node *seChild;
        double xDim[2];
        double yDim[2];
        double centerOfMass;
        double cMassX;
        double cMassY;
        particle *content; //Will be set only and only if the external node
        bool isExternalNode;
        node *parentNode;
	node();
	bool isLeafNode();
	void getNWChild();
	void getNEChild();
	void getSWChild();
	void getSEChild();
};
