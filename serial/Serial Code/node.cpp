#include "node.hpp"


node::node()
{
	nwChild = NULL;
	neChild = NULL;
	swChild = NULL;
	seChild = NULL;
	centerOfMass = 0.0;
	content = NULL;
	isExternalNode = false;
	parentNode = NULL;
}

bool node::isLeafNode()
{
	if (nwChild || neChild || swChild || seChild)
	{
		return false;
	}
	return true;
}

void node::getNWChild()
{
	nwChild = new node();
	nwChild->xDim[0] = xDim[0];
	nwChild->xDim[1] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
	nwChild->yDim[0] = yDim[0] + (yDim[1] - yDim[0]) / 2.0;
	nwChild->yDim[1] = yDim[1];
	nwChild->parentNode = this;
}

void node::getNEChild()
{
	neChild = new node();
	neChild->xDim[0] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
	neChild->xDim[1] = xDim[1];
	neChild->yDim[0] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
	neChild->yDim[1] = yDim[1];
	neChild->parentNode = this;
}

void node::getSWChild()
{
	swChild = new node();
	swChild->xDim[0] = xDim[0];
	swChild->xDim[1] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
	swChild->yDim[0] = yDim[0];
	swChild->yDim[1] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
	swChild->parentNode = this;
}

void node::getSEChild()
{
	seChild = new node();
	seChild->xDim[0] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
	seChild->xDim[1] = xDim[1];
	seChild->yDim[0] = yDim[0];
	seChild->yDim[1] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
	seChild->parentNode = this;
}
