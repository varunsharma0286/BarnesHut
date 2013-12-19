#include "barnesTree.hpp"


barnesTree::barnesTree() //Default must never be invoked
{
	cout<<"\n Default Invoked ";
        m_noOfParticles = 10;
        m_noOfSteps = 1;
        m_theta = THETA;
}

barnesTree::barnesTree(long noOfParticles, long noOfSteps, float theta)
{
	m_noOfParticles = noOfParticles;
	m_noOfSteps = noOfSteps;
	m_theta = theta;
	m_root = new node();
	m_root->xDim[START] = XMIN;
	m_root->xDim[END] = XDIM;
	m_root->yDim[START] = YMIN;
	m_root->yDim[END] = XDIM;
	particleList = new particle*[m_noOfParticles];
	logFile.open("log.txt");
}

void barnesTree::prepareInitList()
{
	srand(time(NULL));
	#pragma omp parallel for
	for (int i = 0; i < m_noOfParticles; i++)
	{
		particleList[i] = new particle();
		particleList[i]->generateXPosn();
		particleList[i]->generateYPosn();
		//cout << "\n" << i << " xPosn = " << particleList[i]->xPosn << " yPosn = " << particleList[i]->yPosn;
	}
}

bool barnesTree::run()
{
	bool largeSet = false;
	prepareInitList();
	for (int i = 0; i < m_noOfParticles; i++)
	{
		insert(m_root,i);
	}

	if(m_noOfParticles > 1000)
	{
		largeSet = true;
		parallelCalcCenterOfMass();
	}
	else
	{
		calcCenterOfMass(m_root);
	}
	logParticlePosition(0);
	for (int steps = 0; steps < m_noOfSteps; steps++)
	{
		#pragma omp parallel for
		for (int idx = 0; idx < m_noOfParticles; idx++)
		{
			particleList[idx]->force = 0.0;
			calcForce(m_root, idx);
			particleList[idx]->updateParticle();
		}
		rebuildTree();
		if(largeSet)
			parallelCalcCenterOfMass();
		else
			calcCenterOfMass(m_root);
		//logParticlePosition(steps+1);
	}
        logParticlePosition(m_noOfSteps);
	return true;
}

void barnesTree::parallelCalcCenterOfMass()
{
	float mass =0.0,posnX=0.0, posnY=0.0;
	#pragma omp parallel sections reduction(+:mass,posnX,posnY)
	{
    		#pragma omp section
    		{	
			if(m_root->nwChild)
			{
				calcCenterOfMass(m_root->nwChild);
				mass = m_root->nwChild->centerOfMass;
				posnX = m_root->nwChild->centerOfMass*m_root->nwChild->cMassX;
				posnY = m_root->nwChild->centerOfMass*m_root->nwChild->cMassY;
			}
    		}
		#pragma omp section
    		{ 	
			if(m_root->neChild)
			{
			        calcCenterOfMass(m_root->neChild);
                                mass = m_root->neChild->centerOfMass;
                                posnX = m_root->neChild->centerOfMass*m_root->neChild->cMassX;
                                posnY = m_root->neChild->centerOfMass*m_root->neChild->cMassY;
			}
    		}
                #pragma omp section
                {
			if(m_root->swChild)
			{
                        	calcCenterOfMass(m_root->swChild);
                                mass = m_root->swChild->centerOfMass;
                                posnX = m_root->swChild->centerOfMass*m_root->swChild->cMassX;
                                posnY = m_root->swChild->centerOfMass*m_root->swChild->cMassY;
			}
                }
                #pragma omp section
                {
			if(m_root->seChild)
			{
	                        calcCenterOfMass(m_root->seChild);
                                mass = m_root->seChild->centerOfMass;
                                posnX = m_root->seChild->centerOfMass*m_root->seChild->cMassX;
                                posnY = m_root->seChild->centerOfMass*m_root->seChild->cMassY;
			}
                }
		
	}
	
	//Calc the center of mass of the root;
	m_root->centerOfMass = mass;
	m_root->cMassX = posnX/mass;
	m_root->cMassY = posnY/mass;
}

void barnesTree::rebuildTree()
{
	clearTree(m_root);
	//cout << "\n Root also got deleted :( Hence creating a new root :P ";
	m_root = new node();
	m_root->xDim[START] = XMIN;
	m_root->xDim[END] = XDIM;
	m_root->yDim[START] = YMIN;
	m_root->yDim[END] = XDIM;
	for (int i = 0; i<m_noOfParticles; i++)
	{
		if(particleList[i]->hasLeftSystem == false)
		{
			insert(m_root, i);
		}
	}
}

void barnesTree::clearTree(node *currNode)
{
	if (currNode)
	{
		if (currNode->isLeafNode())
		{
			delete currNode;
			return;
		}	
	}
	if (currNode->nwChild)
	{
		clearTree(currNode->nwChild);
		currNode->nwChild = NULL;
	}
	if (currNode->neChild)
	{
		clearTree(currNode->neChild);
		currNode->neChild = NULL;
	}
	if (currNode->swChild)
	{
		clearTree(currNode->swChild);
		currNode->swChild = NULL;
	}
	if (currNode->seChild)
	{
		clearTree(currNode->seChild);
		currNode->seChild = NULL;
	}
	delete currNode;
	currNode = NULL;
}

void barnesTree::insert(node *currNode,const int idx)
{
	if (!currNode)
	{
		cout << "\n failed during insertion\n";
		return;
	}
	if (currNode->isLeafNode())
	{
		if (!currNode->content) //Place the particle in this node
		{
			currNode->content = particleList[idx];
			return;
		}
		else if (currNode->content)
		{
			splitNode(currNode, idx);
			currNode->content = NULL;
			insert(currNode, idx);
			return;
		}
	}
	else
	{
		//Check the quadrant to be moved into
		int region = checkQuadrant(currNode, particleList[idx]);
		switch (region)
		{
		case NW:
			if (currNode->nwChild)
			{
				insert(currNode->nwChild, idx);
			}
			else
			{
				//Set it as the nwChild of the currNode
				currNode->getNWChild();
				currNode->nwChild->content = particleList[idx];
			}
			break;
		case NE:
			if (currNode->neChild)
			{
				insert(currNode->neChild, idx);
			}
			else
			{
				//Set it as the neChild of the currNode
				currNode->getNEChild();
				currNode->neChild->content = particleList[idx];
			}
			break;
		case SW:
			if (currNode->swChild)
			{
				insert(currNode->swChild, idx);
			}
			else
			{
				//Set it as the nwChild of the currNode
				currNode->getSWChild();
				currNode->swChild->content = particleList[idx];
			}
			break;
		case SE:
			if (currNode->seChild)
			{
				insert(currNode->seChild, idx);
			}
			else
			{
				//Set it as the nwChild of the currNode
				currNode->getSEChild();
				currNode->seChild->content = particleList[idx];
			}
			break;
		default:
			cout << "\n Entered into the default option: Something went wrong";
			return;
		}
	}
}

void barnesTree::calcCenterOfMass(node *currNode)
{
	if (currNode->isLeafNode())
	{
		if(currNode->content)
		{
			currNode->centerOfMass = MASS;
			currNode->cMassX = currNode->content->xPosn;
			currNode->cMassY = currNode->content->yPosn;
		}
		else
		{
			//cout<<"\n null content";
		}
		return;
	}
	float totalMass = 0.0;
	float cxPosn = 0.0;
	float cyPosn = 0.0;
	if (currNode->nwChild)
	{
		calcCenterOfMass(currNode->nwChild);
		totalMass += currNode->nwChild->centerOfMass;
		cxPosn += (currNode->nwChild->centerOfMass*currNode->nwChild->cMassX);
		cyPosn += (currNode->nwChild->centerOfMass*currNode->nwChild->cMassY);
	}

	if (currNode->neChild)
	{
		calcCenterOfMass(currNode->neChild);
		totalMass += currNode->neChild->centerOfMass;
		cxPosn += (currNode->neChild->centerOfMass*currNode->neChild->cMassX);
		cyPosn += (currNode->neChild->centerOfMass*currNode->neChild->cMassY);
	}

	if (currNode->swChild)
	{
		calcCenterOfMass(currNode->swChild);
		totalMass += currNode->swChild->centerOfMass;
		cxPosn += (currNode->swChild->centerOfMass*currNode->swChild->cMassX);
		cyPosn += (currNode->swChild->centerOfMass*currNode->swChild->cMassY);
	}
	
	if (currNode->seChild)
	{
		calcCenterOfMass(currNode->seChild);
		totalMass += currNode->seChild->centerOfMass;
		cxPosn += (currNode->seChild->centerOfMass*currNode->seChild->cMassX);
		cyPosn += (currNode->seChild->centerOfMass*currNode->seChild->cMassY);
	}
	currNode->centerOfMass = totalMass;
	currNode->cMassX = cxPosn / totalMass;
	currNode->cMassY = cyPosn / totalMass;
}	
void barnesTree::calcForce(node *currNode,long int idx)
{
	//(s/d<Theta)->Where where s is the width of the region represented by the internal node, and d is the 
	//distance between the body and the node's center-of-mass
	if (currNode->isLeafNode())
	{
		if (currNode->content && (currNode->content != particleList[idx]))
		{
			float xDiff = currNode->cMassX - particleList[idx]->xPosn;
			float yDiff = currNode->cMassY - particleList[idx]->yPosn;
			float distance = sqrt((xDiff*xDiff) + (yDiff*yDiff));
			particleList[idx]->forceX += ((G * MASS *MASS) / distance*distance*distance)*xDiff;
                        particleList[idx]->forceY += ((G * MASS *MASS) / distance*distance*distance)*yDiff;
		}
	}
	else
	{
		//Not a leaf node so take the s/d into consideration for this node
		float width = currNode->xDim[1] - currNode->xDim[0];
		float xDiff = currNode->cMassX - particleList[idx]->xPosn;
		float yDiff = currNode->cMassY - particleList[idx]->yPosn;
		float distance = sqrt((xDiff*xDiff) + (yDiff*yDiff));
		if ((width/distance) < THETA)
		{
			particleList[idx]->forceX += (G*MASS*currNode->centerOfMass*xDiff) / (distance*distance*distance);
			particleList[idx]->forceY += (G*MASS*currNode->centerOfMass*yDiff) / (distance*distance*distance);
		}
		else
		{
			if (currNode->nwChild)
			{
				calcForce(currNode->nwChild,idx);
			}
			if (currNode->neChild)
			{
				calcForce(currNode->neChild,idx);
			}
			if (currNode->swChild)
			{
				calcForce(currNode->swChild,idx);
			}
			if (currNode->seChild)
			{
				calcForce(currNode->seChild, idx);
			}
		}
	}
}

bool barnesTree::checkSDRatio(node *currNode, particle *currParticle)
{
	float width = currNode->xDim[1] - currNode->xDim[0];
	float xDiff = currNode->cMassX - currParticle->xPosn;
	float yDiff = currNode->cMassY - currParticle->yPosn;
	float distance = sqrt((xDiff*xDiff) + (yDiff*yDiff));
	if ((distance / width) < THETA)
	{
		return true;
	}
}

int barnesTree::checkQuadrant(node *currNode, particle *temp)
{
	//Check the NW quadrant of the current node
	float xMin = currNode->xDim[0];
	float xMax = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
	float yMin = currNode->yDim[0] + (currNode->yDim[1] - currNode->yDim[0]) / 2.0;
	float yMax = currNode->yDim[1];
	if (temp->xPosn >= xMin && temp->xPosn <= xMax)
	{
		if (temp->yPosn > yMin && temp->yPosn <= yMax)
		{
			return NW;
		}
	}

	//Check the NE Quadrant of the current node
	xMin = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
	xMax = currNode->xDim[1];
	yMin = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);
	yMax = currNode->yDim[1];

	if (temp->xPosn > xMin && temp->xPosn <= xMax)
	{
		if (temp->yPosn > yMin && temp->yPosn <= yMax)
		{
			return NE;
		}
	}

	//Check the SW Quadrant of the current node
	xMin = currNode->xDim[0];
	xMax = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
	yMin = currNode->yDim[0];
	yMax = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);
	if (temp->xPosn >= xMin && temp->xPosn <= xMax)
	{
		if (temp->yPosn >= yMin && temp->yPosn <= yMax)
		{
			return SW;
		}
	}

	//Check the SE Quadrant of the current node
	xMin = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
	xMax = currNode->xDim[1];
	yMin = currNode->yDim[0];
	yMax = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);
	if (temp->xPosn > xMin && temp->xPosn <= xMax)
	{
		if (temp->yPosn >= yMin && temp->yPosn <= yMax)
		{
			return SE;
		}
	}

	//Must not reach till this point
	//If here that means an error
	return INVALID;
}

void barnesTree::splitNode(node *currNode, int idx)
{
	//Here we can assume that the currNode is the leaf node but already has the content
	//So we need to move both the contents down
	if (!currNode)
	{
		cout << "\n Error::splitNode() got currNode as null-> returning without doing anythig";
		return;
	}
	if (currNode->xDim[1] - currNode->xDim[0] == 1.0f)
	{
		//Cannot split this node further.
		//cout << "\n Cannot split this node further.";
		currNode->centerOfMass += particleList[idx]->mass;
		return;
	}
	int contentQuad = checkQuadrant(currNode, currNode->content);
	switch (contentQuad)
	{
	case NW:
		currNode->getNWChild();
		currNode->nwChild->content = currNode->content;
		break;
	case NE:
		currNode->getNEChild();
		currNode->neChild->content = currNode->content;
		break;
	case SW:
		currNode->getSWChild();
		currNode->swChild->content = currNode->content;
		break;
	case SE:
		currNode->getSEChild();
		currNode->seChild->content = currNode->content;
		break;
	}
}

void barnesTree::logParticlePosition(int stepNumber)
{
	int i =  m_noOfParticles/4;
	stringstream output[4];
	output[0].precision(8);
	output[1].precision(8);
	output[2].precision(8);
	output[3].precision(8);

	#pragma omp parallel sections
	{
		#pragma omp section
		for(int j = 0;j<i;j++)
		{
			output[0]<<"\n" << j << " xPosn = "<< particleList[j]->xPosn << " yPosn = " <<particleList[j]->yPosn;
			output[0]<<" has left= "<<particleList[j]->hasLeftSystem;
		}
		#pragma omp section
                for(int j = i;j<2*i;j++)
                {
                        output[1]<<"\n" << j<< " xPosn = " << particleList[j]->xPosn << " yPosn = " << particleList[j]->yPosn;
			output[1]<<" has left= "<<particleList[j]->hasLeftSystem;
                }
                #pragma omp section
                for(int j =2*i;j<3*i;j++)
                {
                        output[2]<<"\n" << j << " xPosn = " << particleList[j]->xPosn << " yPosn = " << particleList[j]->yPosn;
			output[2]<<" has left= "<<particleList[j]->hasLeftSystem;
                }
                #pragma omp section
                for(int j =3*i;j<m_noOfParticles;j++)
                {
                        output[3]<<"\n" << j << " xPosn = " << particleList[j]->xPosn << " yPosn = " << particleList[j]->yPosn;
			output[3]<<" has left="<<particleList[j]->hasLeftSystem;
                }
	}

	stringstream fileName;
	int noOfThreads = 0;
        #pragma omp parallel
        {
                #pragma omp master
                {
                        #ifdef _OPENMP
                        noOfThreads = omp_get_num_threads();
                        #endif
                }
        }

	fileName<<"output/step_"<<noOfThreads<<"_"<<m_noOfParticles<<"_"<<stepNumber<<".txt";
	string title = fileName.str();
	//cout<<"\n fileName = "<<fileName.str();
	ofstream test(title.c_str());
	if(!test)
		cout<<"\n error creating the output file";
	else
	{
		test<<output[0].str()<<output[1].str()<<output[2].str()<<output[3].str();
		test.close();
	}
}

