#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args)
{
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
    gamma = params.k;
    Rc = params.r_cut;
    ttotal = params.ttime;
    nsteps = (int) ttotal / dt;
    
    ncells = params.n_particles;
    setIntegrator(params.integrator_a);
    
    cout << "dt: " << dt << endl;
    cout << "a: "  << a <<  endl;
    cout << "d: " << d << endl;
}

Simulator::Simulator(const Simulator& orig) 
{
}

Simulator::~Simulator() 
{
}

void Simulator::addCell(const Cell& newCell)
{
    cells.push_back(newCell);
    numberofCells++;
}

void Simulator::addCell()
{
   SimpleTriangulation sm(params.d);
   list<Triangle> tris = sm.triangulate();
   Cell newCell(tris);
   addCell(newCell);
}

void Simulator::calcForces()
{
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces();
    }
    
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < numberofCells; j++) 
        {
            if (i != j )
                cells[i].calcForces(cells[j]);
        }
    }
}

void Simulator::integrateEuler()
{
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberV; j++)
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
        for (int j = 0; j < cells[i].numberV; j++)
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
        for (int j = 0; j < cells[i].numberV; j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
            calcForces();
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; // v(t+1) = v(t+1)_1 + 0.5 * dt * a(t+1)
        }
    
    }    
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

void Simulator::saveCellsState()
{
    int index;
    cout << "saving in: " << params.output_file << endl;
    ofstream os(params.output_file);
    for (int i = 0; i < numberofCells; i++)
    {
        os << cells[i].numberV << "\n" ;
        for (int j = 0; j < cells[i].numberV; j++)
        {
            index = (cells[i].vertices[j].getId()+1) ;
            os << "H" << index << " "<< cells[i].vertices[j].xyz.x << " " << cells[i].vertices[j].xyz.y << " " << cells[i].vertices[j].xyz.z << "\n";
        }  
    }
    os.close();
}
