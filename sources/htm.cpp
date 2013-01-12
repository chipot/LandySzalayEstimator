#include "../includes/htm.hpp"

HTM* ICoDF_HTM::HTM::_singleton = NULL; // Initialize the singleton to NULL.

// ADD POINT
void ICoDF_HTM::HTM::AddPoint(const double& ra, const double& dec)
{
  PointInfo_t* info = new PointInfo_t;
  info->_ra = ra;
  info->_dec = dec;
  if (this->SelectRootTrixel(info))
    this->_pointList.push(info);
  else
    delete info;
}

// CreateHTM
bool ICoDF_HTM::HTM::CreateHTM()
{
  PointInfo_t* pt;
  while (!this->_pointList.empty())
    {
      pt = _pointList.front();
      this->_pointList.pop();
      this->AssignPoint(pt);
    }
  return true;
}

// AssignPoint
std::string ICoDF_HTM::HTM::AssignPoint(PointInfo_t* pt)
{
  if (pt->_current->_nbChildObject > 1)
    {
      unsigned short int index = GetIndex(pt->_current, pt);
      if (index == (unsigned int short)~0)
	{
	  return std::string("error");
	}

      if (pt->_current->_children[index] != NULL)
	{
  pt->_current = pt->_current->_children[index];

	  this->_pointList.push(pt);
	}
      else
	{
	  pt->_current->_children[index] = CreateTrixelChild(pt->_current, index);
	  pt->_current = pt->_current->_children[index];
	  this->_pointList.push(pt);
	}
    }
  else if (pt->_current->_nbChildObject == 1)
    {
      ++pt->_current->_nbChildObject;
      unsigned short int indexCurrent = GetIndex(pt->_current, pt);
      PointInfo_t* old = pt->_current->_info;
      pt->_current->_info = NULL;
      unsigned short int indexOld = GetIndex(old->_current, old);
      if (indexCurrent == (unsigned short int)~0)
	{
	  return std::string("");
	}
      if (indexOld == (unsigned short int)~0)
	{
	  return std::string("");
	}
      if (pt->_current->_children == NULL)
	{
	  CreateTrixelChildren(pt->_current);
	  CreateTrixelChild(pt->_current, indexOld);
	  if (indexOld != indexCurrent)
	    CreateTrixelChild(pt->_current, indexCurrent);
	  auto it = this->_points.find(pt->_current->_HTMId);
	  this->_points.erase(it);
	}
      old->_current = pt->_current->_children[indexOld];
      this->_pointList.push(pt);
      this->_pointList.push(old);
      pt->_current = pt->_current->_children[indexCurrent];
    }
  else if (pt->_current->_nbChildObject == 0)
    {
      pt->_current->_info = pt;
      this->_points[pt->_current->_HTMId] = pt;
      pt->_current->_nbChildObject = 1;
      return pt->_current->_HTMId;
    }
  return pt->_current->_HTMId;
}

double ICoDF_HTM::HTM::getMinRa(void)
{
  return this->_raQueue.top();
}

double ICoDF_HTM::HTM::getMinDec(void)
{
  return this->_decQueue.top();
}

double ICoDF_HTM::HTM::getMaxRa(void)
{
  double raMax;
  while (!this->_raQueue.empty())
    {
      raMax = this->_raQueue.top();
      this->_raQueue.pop();
    }
  return raMax;
}

double ICoDF_HTM::HTM::getMaxDec(void)
{
  double decMax;
  while (!this->_decQueue.empty())
    {
      decMax = this->_decQueue.top();
      this->_decQueue.pop();
    }
  return decMax;
}

void ICoDF_HTM::HTM::itemsToStore(const double& ra, const double& dec)
{
  this->_raQueue.push(ra);
  this->_decQueue.push(dec);
}

HTMConstraint_t* ICoDF_HTM::HTM::SetConstraint(PointInfo_t* pt, double& radius)
{
  HTMConstraint_t* constraint = NULL;
  if (pt != NULL)
    {
      
      if (IsCorrectRA(pt->_ra) && IsCorrectDEC(pt->_dec))
	{
	  double rProjection = sin(90 - abs(pt->_dec));
	  double x = rProjection * cos(pt->_ra);
	  double y = rProjection * sin(pt->_ra);
	  double z = cos(90 - abs(pt->_dec));
	  Eigen::Vector3d p(x, y, z);
	  constraint = new HTMConstraint_t;
	  static std::queue<trixel_t*> workingList;
	  std::queue<trixel_t*>().swap(workingList);

	  for(unsigned int i = 0; i < 4; ++i)
	    {
	      unsigned short int inside = 0;
	      if (this->_octahedron->_rootTrixels[i] != NULL)
		{
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[0]) > radius)
		    ++inside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[1]) > radius)
		    ++inside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[2]) > radius)
		    inside++;
		}
	      if (inside == 3)
		constraint->_inside.push_back(this->_octahedron->_rootTrixels[i]);
	      else if (inside > 0)
		workingList.push(this->_octahedron->_rootTrixels[i]);
	      else
		{
		  Eigen::Vector3d tmpVec1 = _octahedron->_rootTrixels[i]->_vertices[1] - _octahedron->_rootTrixels[i]->_vertices[0];
		  Eigen::Vector3d tmpVec2 = _octahedron->_rootTrixels[i]->_vertices[2] - _octahedron->_rootTrixels[i]->_vertices[1];
		  Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2);
		  Eigen::Vector3d trixelBoundary = tmpVec3 / tmpVec3.norm();

		  double theta = acos(trixelBoundary.dot(p) / (trixelBoundary.norm() * p.norm()));
		  double phi1 = acos(trixelBoundary.dot(Eigen::Vector3d(1,0,0)) / (trixelBoundary.norm()));
		  double phi2 = acos(p.dot(Eigen::Vector3d(1,0,0)) / p.norm());
		  if (theta < phi1 + phi2)
		    {
		      if (!(_octahedron->_rootTrixels[i]->_vertices[0].cross(_octahedron->_rootTrixels[i]->_vertices[1]).dot(p) < 0 &&
			    _octahedron->_rootTrixels[i]->_vertices[1].cross(_octahedron->_rootTrixels[i]->_vertices[2]).dot(p) < 0 &&
			    _octahedron->_rootTrixels[i]->_vertices[2].cross(_octahedron->_rootTrixels[i]->_vertices[0]).dot(p)))
			{
			  constraint->_partial.push_back(_octahedron->_rootTrixels[i]);
			}
		    }
		}      
	    }
	  while (workingList.size() > 0)
	    {
	      trixel_t* tmp = workingList.front();
	      workingList.pop();
	      unsigned short int inside = 0;
	      if (p.dot(tmp->_vertices[0]) > radius)
		++inside;
	      if (p.dot(tmp->_vertices[2]) > radius)
		++inside;
	      if (p.dot(tmp->_vertices[1]) > radius)
		++inside;
	      
	      if (inside == 3)
		constraint->_inside.push_back(tmp);
	      else if (inside > 0)
		{
		  if (tmp->_children != NULL)
		    {
		      for (unsigned int i = 0; i < 4; ++i)
			if (tmp->_children[i] != NULL)
			  workingList.push(tmp->_children[i]);
		    }
		  else
		    constraint->_partial.push_back(tmp);
		}
	      else
		{
		  Eigen::Vector3d tmpVec1 = tmp->_vertices[1] - tmp->_vertices[0];
		  Eigen::Vector3d tmpVec2 = tmp->_vertices[2] - tmp->_vertices[1];
		  Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2);
		  Eigen::Vector3d trixelBoundary = tmpVec3 / tmpVec3.norm();

		  double theta = acos(trixelBoundary.dot(p) / (trixelBoundary.norm() * p.norm()));
		  double phi1 = acos(trixelBoundary.dot(Eigen::Vector3d(1,0,0)) / (trixelBoundary.norm()));
		  double phi2 = acos(p.dot(Eigen::Vector3d(1,0,0)) / p.norm());
		  if (theta < phi1 + phi2)
		    {
		      if (!(tmp->_vertices[0].cross(tmp->_vertices[1]).dot(p) < 0 &&
			    tmp->_vertices[1].cross(tmp->_vertices[2]).dot(p) < 0 &&
			    tmp->_vertices[2].cross(tmp->_vertices[0]).dot(p)))
			{
			  constraint->_partial.push_back(tmp);
			}
		    }
		}
	    } 
	}
      else
	{
	  LS_ADDMSG(LogService::WARNING, "ICoDF_HTM::HTM::SetConstraint", "Given right ascension and/or declination has incorrect value");
	  return NULL;
	}
    }
  else
    {
      LS_ADDMSG(LogService::NOTICE, "HTM::SetConstraint", "Given <pt> has a NULL value");
    }

  return constraint;
}

inline std::pair<double, double>       ICoDF_HTM::HTM::CalcCoordPoint(std::pair<double, double>& a, std::pair<double, double>& b)
{
  std::pair<double, double>	result;

  result.first = a.first - b.first;
  result.second = a.second - b.second;
  return result;
}

bool ICoDF_HTM::HTM::SelectRootTrixel(PointInfo_t* pt)
{
  for (int i = 0; i < 8; ++i)
    {
      Eigen::Vector3d* v = this->_octahedron->_rootTrixels[i]->_vertices;
      double rProjection = sin(90 - abs(pt->_dec));
      double x = rProjection * cos(pt->_ra);
      double y = rProjection * sin(pt->_ra);
      double z = cos(90 - abs(pt->_dec));
      Eigen::Vector3d p(x, y, z);

      if (v[0].cross(v[1]).dot(p) > 0 &&
	  v[1].cross(v[2]).dot(p) > 0 &&
	  v[2].cross(v[0]).dot(p) > 0)
	{
	  pt->_current = this->_octahedron->_rootTrixels[i];
	  return true;
	}
    }
  return false;
}

inline double				ICoDF_HTM::HTM::Scal(std::pair<double, double>& v1, std::pair<double, double>& v2) const
{
  return ((v1.first * v2.first) + (v1.second * v2.second));
}

bool				ICoDF_HTM::HTM::CheckPointInTriangle(std::pair<double, double> A,
						     std::pair<double, double> B,
						     std::pair<double, double> C,
						     std::pair<double, double> P)
{
  std::pair<double, double> 	v0 = CalcCoordPoint(C, A);
  std::pair<double, double> 	v1 = CalcCoordPoint(B, A);
  std::pair<double, double> 	v2 = CalcCoordPoint(P, A);

  const double			dot00 = Scal(v0, v0);
  const double			dot01 = Scal(v0, v1);
  const double			dot02 = Scal(v0, v2);
  const double			dot11 = Scal(v1, v1);
  const double			dot12 = Scal(v1, v2);

  const double       		invDenom = 1 / (dot00 * dot11 - dot01 * dot01);

  const double       		u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  const double       		v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  return ((u >= 0) && (v >= 0) && (u + v < 1));
}

/// Create the instance of the HTM core if applicable and/or return a pointer to it.
HTM* ICoDF_HTM::HTM::GetInstance()
{
  if (HTM::_singleton == NULL)
    {
      HTM::_singleton = new HTM();
    }
  return HTM::_singleton;
}

/// DELETE
/// Delete the instance of the HTM core
void ICoDF_HTM::HTM::Delete()
{
  if (HTM::_singleton != NULL)
    {
      delete HTM::_singleton;
      HTM::_singleton = NULL;
    }
}

/// TwoPointsCorrelation
unsigned int ICoDF_HTM::HTM::TwoPointsCorrelation(double& radius, double& delta)
{
  unsigned int nbPairs = 0;
  std::map<std::string, PointInfo_t*>::iterator it;

  double infLimit = radius - delta;
  if (infLimit < 0) infLimit = 0;
  double supLimit = radius + delta;
  HTMConstraint_t *constraint = new HTMConstraint_t;

  for (it = this->_points.begin(); it != this->_points.end(); ++it)
    {
      PointInfo_t* pt = (*it).second;
      if (IsCorrectRA(pt->_ra) && IsCorrectDEC(pt->_dec))
	{
	  double rProjection = sin(90 - abs(pt->_dec));
	  double x = rProjection * cos(pt->_ra);
	  double y = rProjection * sin(pt->_ra);
	  double z = cos(90 - abs(pt->_dec));
	  Eigen::Vector3d p(x, y, z);

	  static std::queue<trixel_t*> workingList;
	  std::queue<trixel_t*>().swap(workingList);

	  for (unsigned int i = 0; i < 4; ++i)
	    {
	      unsigned short int infInside = 0;
	      unsigned short int supInside = 0;
	      if (this->_octahedron->_rootTrixels[i] != NULL)
		{
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[0]) > infLimit)
		    ++infInside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[1]) > infLimit)
		    ++infInside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[2]) > infLimit)
		    infInside++;

		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[0]) > supLimit)
		    ++supInside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[1]) > supLimit)
		    ++supInside;
		  if (p.dot(this->_octahedron->_rootTrixels[i]->_vertices[2]) > supLimit)
		    supInside++;
		}
	      if (supInside == 3 && infInside == 0)
		constraint->_inside.push_back(this->_octahedron->_rootTrixels[i]);
	      if ((supInside == 3 && infInside > 0) || supInside > 0)
		workingList.push(this->_octahedron->_rootTrixels[i]);
	      else
		{
		  Eigen::Vector3d tmpVec1 = _octahedron->_rootTrixels[i]->_vertices[1] - _octahedron->_rootTrixels[i]->_vertices[0];
		  Eigen::Vector3d tmpVec2 = _octahedron->_rootTrixels[i]->_vertices[2] - _octahedron->_rootTrixels[i]->_vertices[1];
		  Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2);
		  Eigen::Vector3d trixelBoundary = tmpVec3 / tmpVec3.norm();

		  double theta = acos(trixelBoundary.dot(p) / (trixelBoundary.norm() * p.norm()));
		  double phi1 = acos(trixelBoundary.dot(Eigen::Vector3d(1,0,0)) / (trixelBoundary.norm()));
		  double phi2 = acos(p.dot(Eigen::Vector3d(1,0,0)) / p.norm());
		  if (theta < phi1 + phi2)
		    {
		      if (!(_octahedron->_rootTrixels[i]->_vertices[0].cross(_octahedron->_rootTrixels[i]->_vertices[1]).dot(p) < 0 &&
			    _octahedron->_rootTrixels[i]->_vertices[1].cross(_octahedron->_rootTrixels[i]->_vertices[2]).dot(p) < 0 &&
			    _octahedron->_rootTrixels[i]->_vertices[2].cross(_octahedron->_rootTrixels[i]->_vertices[0]).dot(p)))
			{
			  constraint->_partial.push_back(_octahedron->_rootTrixels[i]);
			}
		    }		  
		}
	    }
	  while (workingList.size() > 0)
	    {
	      trixel_t* tmp = workingList.front();
	      workingList.pop();
	      unsigned short int infInside = 0;
	      unsigned short int supInside = 0;
	      if (p.dot(tmp->_vertices[0]) > infLimit)
		++infInside;
	      if (p.dot(tmp->_vertices[2]) > infLimit)
		++infInside;
	      if (p.dot(tmp->_vertices[1]) > infLimit)
		++infInside;

	      if (p.dot(tmp->_vertices[0]) > supLimit)
		++supInside;
	      if (p.dot(tmp->_vertices[2]) > supLimit)
		++supInside;
	      if (p.dot(tmp->_vertices[1]) > supLimit)
		++supInside;
	      
	      if (supInside == 3 && infInside == 0)
		constraint->_inside.push_back(tmp);
	      else if ((supInside == 3 && infInside > 0) 
		       || supInside > 0)
		{
		  if (tmp->_children != NULL)
		    {
		      for (unsigned int i = 0; i < 4; ++i)
			if (tmp->_children[i] != NULL)
			  workingList.push(tmp->_children[i]);
		    }
		  else
		    constraint->_partial.push_back(tmp);
		}
	      else
		{
		  Eigen::Vector3d tmpVec1 = tmp->_vertices[1] - tmp->_vertices[0];
		  Eigen::Vector3d tmpVec2 = tmp->_vertices[2] - tmp->_vertices[1];
		  Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2);
		  Eigen::Vector3d trixelBoundary = tmpVec3 / tmpVec3.norm();
		  
		  double theta = acos(trixelBoundary.dot(p) / (trixelBoundary.norm() * p.norm()));
		  double phi1 = acos(trixelBoundary.dot(Eigen::Vector3d(1,0,0)) / (trixelBoundary.norm()));
		  double phi2 = acos(p.dot(Eigen::Vector3d(1,0,0)) / p.norm());
		  if (theta < phi1 + phi2)
		    {
		      if (!(tmp->_vertices[0].cross(tmp->_vertices[1]).dot(p) < 0 &&
			    tmp->_vertices[1].cross(tmp->_vertices[2]).dot(p) < 0 &&
			    tmp->_vertices[2].cross(tmp->_vertices[0]).dot(p)))
			{
			  constraint->_partial.push_back(tmp);
			}
		    }
		}
	    }
	}
    }

  for (auto it2 = constraint->_inside.begin(); it2 != constraint->_inside.end(); ++it2)
    {
      nbPairs += (*it2)->_nbChildObject;
    }
  for (auto it2 = constraint->_partial.begin(); it2 != constraint->_partial.end(); ++it2)
    {
      nbPairs += 1;
    }

  delete constraint;
  return nbPairs;
}

/// SelectOctahedronTrixel
trixel_t* ICoDF_HTM::HTM::SelectRootOctahedronTrixel(const double& ra, const double& dec)
{
  unsigned int trixelNumber = 0;

  if (CheckPointInTriangle(std::make_pair(0.0, 0.0), std::make_pair(90.0, 0.0), std::make_pair(45.0, 90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 0;
    }
  else if (CheckPointInTriangle(std::make_pair(90.0, 0.0), std::make_pair(180.0, 0.0), std::make_pair(135.0, 90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 1;
    }
  else if (CheckPointInTriangle(std::make_pair(180.0, 0.0), std::make_pair(270.0, 0.0), std::make_pair(225.0, 90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 2;
    }
  else if (CheckPointInTriangle(std::make_pair(270.0, 0.0), std::make_pair(360.0, 0.0), std::make_pair(315.0, 90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 3;
    }
  else if (CheckPointInTriangle(std::make_pair(0.0, 0.0), std::make_pair(90.0, 0.0), std::make_pair(45.0, -90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 4;
    }
  else if (CheckPointInTriangle(std::make_pair(90.0, 0.0), std::make_pair(180.0, 0.0), std::make_pair(135.0, -90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 5;
    }
  else if (CheckPointInTriangle(std::make_pair(180.0, 0.0), std::make_pair(270.0, 0.0), std::make_pair(225.0, -90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 6;
    }
  else if (CheckPointInTriangle(std::make_pair(270.0, 0.0), std::make_pair(360.0, 0.0), std::make_pair(315.0, -90.0), std::make_pair(ra, dec)))
    {
      trixelNumber = 7;
    }
  else
    {
      trixelNumber = 8;
      return NULL;
    }
  return this->_octahedron->_rootTrixels[trixelNumber];
}

void	ICoDF_HTM::HTM::CreateOctahedron(void)
{
  this->_octahedron = new Octahedron_t;
  this->_octahedron->_rootTrixels = new trixel_t*[8];
  Eigen::Vector3d v0( 0,  0,  1);
  Eigen::Vector3d v1( 1,  0,  0);
  Eigen::Vector3d v2( 0,  1,  1);
  Eigen::Vector3d v3(-1,  0,  0);
  Eigen::Vector3d v4( 0, -1,  0);
  Eigen::Vector3d v5( 0,  0, -1);

  this->_octahedron->_rootTrixels[0] = CreateRootTrixel(std::string("S0"));
  this->_octahedron->_rootTrixels[0]->_vertices[0] = v1;
  this->_octahedron->_rootTrixels[0]->_vertices[1] = v5;
  this->_octahedron->_rootTrixels[0]->_vertices[2] = v2;
  this->_octahedron->_rootTrixels[1] = CreateRootTrixel(std::string("S1"));
  this->_octahedron->_rootTrixels[1]->_vertices[0] = v2;
  this->_octahedron->_rootTrixels[1]->_vertices[1] = v5;
  this->_octahedron->_rootTrixels[1]->_vertices[2] = v3;
  this->_octahedron->_rootTrixels[2] = CreateRootTrixel(std::string("S2"));
  this->_octahedron->_rootTrixels[2]->_vertices[0] = v3;
  this->_octahedron->_rootTrixels[2]->_vertices[1] = v5;
  this->_octahedron->_rootTrixels[2]->_vertices[2] = v4;
  this->_octahedron->_rootTrixels[3] = CreateRootTrixel(std::string("S3"));
  this->_octahedron->_rootTrixels[3]->_vertices[0] = v4;
  this->_octahedron->_rootTrixels[3]->_vertices[1] = v5;
  this->_octahedron->_rootTrixels[3]->_vertices[2] = v1;
  this->_octahedron->_rootTrixels[4] = CreateRootTrixel(std::string("N0"));
  this->_octahedron->_rootTrixels[4]->_vertices[0] = v1;
  this->_octahedron->_rootTrixels[4]->_vertices[1] = v0;
  this->_octahedron->_rootTrixels[4]->_vertices[2] = v4;
  this->_octahedron->_rootTrixels[5] = CreateRootTrixel(std::string("N1"));
  this->_octahedron->_rootTrixels[5]->_vertices[0] = v4;
  this->_octahedron->_rootTrixels[5]->_vertices[1] = v0;
  this->_octahedron->_rootTrixels[5]->_vertices[2] = v3;
  this->_octahedron->_rootTrixels[6] = CreateRootTrixel(std::string("N2"));
  this->_octahedron->_rootTrixels[6]->_vertices[0] = v3;
  this->_octahedron->_rootTrixels[6]->_vertices[1] = v0;
  this->_octahedron->_rootTrixels[6]->_vertices[2] = v2;
  this->_octahedron->_rootTrixels[7] = CreateRootTrixel(std::string("N3"));
  this->_octahedron->_rootTrixels[7]->_vertices[0] = v2;
  this->_octahedron->_rootTrixels[7]->_vertices[1] = v0;
  this->_octahedron->_rootTrixels[7]->_vertices[2] = v1;
}

void	ICoDF_HTM::HTM::Display(trixel_t* current, std::ofstream& fstream)
{
  if (current != NULL)
    {
      if (current->_children != NULL)
	{
	  for (auto i = 0; i < 4; ++i)
	    {
	      if (current->_children[i] != NULL)
		{
		  Display(current->_children[i], fstream);
		  if (current->_children[0])
		    {
		      if (this->_points[current->_children[0]->_HTMId])
			{
			  PointInfo_t* info = this->_points[current->_children[0]->_HTMId];
			  fstream << "Item stored at trixel : " << current->_children[0]->_HTMId << " with right ascension and declinaison at " << info->_ra << " " << info->_dec << std::endl;
			}
		    }
		  else if (current->_children[1])
		    {
		      if (this->_points[current->_children[1]->_HTMId])
			{
			  PointInfo_t* info = this->_points[current->_children[1]->_HTMId];
			  fstream << "Item stored at trixel : " << current->_children[1]->_HTMId << " with right ascension and declinaison at " << info->_ra << " " << info->_dec << std::endl;
			}
		    }
		  else if (current->_children[2])
		    {
		      if (this->_points[current->_children[2]->_HTMId])
			{
			  PointInfo_t* info = this->_points[current->_children[2]->_HTMId];
			  fstream << "Item stored at trixel : " << current->_children[2]->_HTMId << " with right ascension and declinaison at " << info->_ra << " " << info->_dec << std::endl;
			}
		    }
		  else if (current->_children[3])
		    {
		      if (this->_points[current->_children[3]->_HTMId])
			{
			  PointInfo_t* info = this->_points[current->_children[3]->_HTMId];
			  fstream << "Item stored at trixel : " << current->_children[3]->_HTMId << " with right ascension and declinaison at " << info->_ra << " " << info->_dec << std::endl;
			}
		    }
		}
	    }
	}
    }
}

void	ICoDF_HTM::HTM::FreeAllTrixels(trixel_t* current)
{
  if (current != NULL)
    {
      if (current->_children != NULL)
	{
	  for (auto i = 0; i < 4; ++i)
	    {
	      if (current->_children[i] != NULL)
		{
		  FreeAllTrixels(current->_children[i]);
		  delete current->_children[i];
		}
	    }
	  delete[] current->_children;
	}
    }
}
  
void	ICoDF_HTM::HTM::DeleteOctahedron(void)
{
  std::ofstream fstream;
  fstream.open("log");
  for (auto i = 0; i < 8; ++i)
    this->Display(this->_octahedron->_rootTrixels[i], fstream);
  fstream.close();
  for (auto i = 0; i < 8; ++i)
    this->FreeAllTrixels(this->_octahedron->_rootTrixels[i]);
  for (auto it = this->_points.begin(); it != this->_points.end(); ++it)
    delete it->second;
  for (auto i = 0; i < 8; ++i)
    delete this->_octahedron->_rootTrixels[i];
  this->_points.clear();
  delete this->_octahedron;
}

/// Create the HTM
ICoDF_HTM::HTM::HTM()
{
  LS_ADDMSG(LogService::NOTICE, "HTM", "HTM core created");
}

ICoDF_HTM::HTM::~HTM()
{
  LS_ADDMSG(LogService::NOTICE, "HTM", "HTM core deleted");
}
