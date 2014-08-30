#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args), numberofCells(0), boxSize(DBL_MAX)
{
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
    gamma = params.k;
    Rc = params.r_cut;
    ttotal = params.ttime;
    trajfile = params.traj_file;
    script = params.output_file;
    nsteps = (int) ttotal / dt;
    setIntegrator(params.integrator_a);
}

Simulator::Simulator(const Simulator& orig) 
{
}

Simulator::~Simulator() 
{
}

void Simulator::addCell(const Cell& newCell)
{
    try
    {
        if (cells.size() == MAX_CELLS) 
            throw MaxSizeException();
        
        cells.push_back(newCell);
        numberofCells++;        
    } 
    catch (MaxSizeException& e)
    {
        e.what();
        return;
    }
}

void Simulator::addCell()
{
    
    try
    {
        if (cells.size() == MAX_CELLS) 
            throw MaxSizeException();
        
        SimpleTriangulation sm(params.d);
        list<Triangle> tris = sm.triangulate();
        Cell newCell(tris);
   
        newCell.setA(a);
        newCell.setDp(dp);
        newCell.setRc(Rc);
        newCell.setGamma(gamma);
   
        addCell(newCell);       
    } 
    catch (MaxSizeException& e)
    {
        e.what();
        return;
    }
}

void Simulator::calcForces()
{
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces();
    }
    
    for (int i = 0; i < numberofCells; i++) 
    {
        for (int j = 0; j < numberofCells; j++) 
        {
            if (i != j)
            {
                cells[i].calcForces(cells[j]);    
            }
        }
    }
}

void Simulator::integrateEuler()
{
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].xyz += cells[i].vertices[j].velocity * dt;
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].velocity += cells[i].vertices[j].force * dt / m;
        }
    
    }
}

void Simulator::integrateDampedEuler()
{
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].xyz += cells[i].vertices[j].velocity * dt;
            cells[i].vertices[j].velocity *= 0.0;
        }
    
    }    
}

void Simulator::integrateVv()
{
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
            //calcForces();
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; // v(t+1) = v(t+1)_1 + 0.5 * dt * a(t+1)
        }
    
    }    
}

void Simulator::simulate()
{
    for (int i = 0; i < nsteps; i++)
    {
        doStep();
    }
        
}

void Simulator::simulate(int steps)
{
    renderScript();
    int lastCellIndex = 0;
    int index;
    ofstream os;
    os.open(trajfile);
    os << getTotalVertices() << "\n";
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            index = cells[i].vertices[j].getId()+1 + lastCellIndex;
            os << "H" << index << " ";
            os << cells[i].vertices[j].xyz.x << " " << cells[i].vertices[j].xyz.y << " " << cells[i].vertices[j].xyz.z;
            os << "\n";
        }
        lastCellIndex = cells[i].numberofVertices();
    }

    for (int i = 0; i < steps; i++)
    {
        doStep();
        os << getTotalVertices() << "\n";
        lastCellIndex = 0;
        for (int i = 0; i < numberofCells; i++)
        {
            for (int j = 0; j < cells[i].numberofVertices(); j++)
            {
                index = cells[i].vertices[j].getId()+1 + + lastCellIndex;
                os << "H" << index << " ";
                os << cells[i].vertices[j].xyz.x << " " << cells[i].vertices[j].xyz.y << " " << cells[i].vertices[j].xyz.z;
                os << "\n";
            }
            lastCellIndex = cells[i].numberofVertices();
        } 
    }
    
    os.close();
}

void Simulator::doStep()
{
    calcForces();
    integrate();
}

void Simulator::setIntegrator(void (Simulator::*functoall)())
{
    integrator = functoall;
}

void Simulator::setIntegrator(char* token)
{
    if (STRCMP (token, "vv"))
    {
        this->setIntegrator(&Simulator::integrateVv);
    }
    else if (STRCMP (token, "eu"))
    {
        this->setIntegrator(&Simulator::integrateEuler);
    }
    else if (STRCMP (token, 'de'))
    {
        this->setIntegrator(&Simulator::integrateDampedEuler);
    }
    else
    {
        this->setIntegrator(&Simulator::integrateEuler);
    }
}


void Simulator::integrate()
{
    (*this.*integrator)();
}

void Simulator::moveCell(const Vector3D& v3d, int cellid)
{
    cells[cellid].addXYZ(v3d);
}

void Simulator::addCellVel(const Vector3D& v3d, int cellid)
{
    cells[cellid].addVelocity(v3d);
}

void Simulator::renderScript()
{
    ofstream os;
    os.open(script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "cmd.do(\"load " << trajfile << ", cells\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter elem h, vdw=0.1\")\n";
    os << "cmd.do(\"rebuild\")\n";

    int iidx, jidx;
    int lastCellIndex = 0;
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            iidx = cells[i].vertices[j].getId() + 1 + lastCellIndex;
            for (int k = 0; k < cells[i].vertices[j].nneigh; k++)
            {
                jidx = cells[i].vertices[j].neighbors[k] + 1 + lastCellIndex;
                os << "cmd.do(\"bond /cells///UNK`/H"<< iidx << ", /cells///UNK`/H" << jidx << "\")\n";
            }
            
        }
        lastCellIndex = cells[i].numberofVertices();
    }
    
    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n";
    os.close();
}

void Simulator::saveCellsState(const char* fileout)
{
    int index;
    ofstream os;
    os.open(fileout);
    cout << fileout << endl;
    int totalnumber = getTotalVertices();
   
    os << totalnumber << "\n" ;
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            index = (cells[i].vertices[j].getId()+1) ;
            os << "H" << index << " "<< cells[i].vertices[j].xyz.x << " " << cells[i].vertices[j].xyz.y << " " << cells[i].vertices[j].xyz.z << "\n";
        }  
    }
    os.close();    
}

void Simulator::saveCellsState()
{
    saveCellsState(params.output_file);
}

void Simulator::setBoxSize(double bs)
{
    boxSize = bs;
}

void Simulator::printCell(int i)
{
    this->cells[i].printCell();
}

int Simulator::getTotalVertices()
{
    int totalnumber = 0;
    for (int i = 0; i < numberofCells; i++)
    {
        totalnumber += cells[i].numberofVertices();
    }
    return totalnumber;
}