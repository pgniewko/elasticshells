#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args), numberofCells(0), boxSize(DBL_MAX)
{
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
    bs = params.bs;
    drawBox = params.draw_box;
    pbc = params.pbc;
    gamma = params.k;
    Rc = params.r_cut;
    ttotal = params.ttime;
    trajfile = params.traj_file;
    script = params.output_file;
    nsteps = (int) ttotal / dt;
    setIntegrator(params.integrator_a);
    
    try
    {
        diagnoseParams();
    }
    catch(NotImplementedException& e)
    {
        cout <<  e.what() << endl;
        exit(1);
    }
    catch(DataException& e)
    {
        cout << e.what() << endl;
        exit(1);
    }
}

Simulator::Simulator(const Simulator& orig) {}

Simulator::~Simulator() {}

void Simulator::diagnoseParams()
{
    if (d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                "Single point representation is not implemented yet. "
                "Simulator is about to terminate !");
    if (d > 7)
        throw DataException("DataException:\n"
                "Depth of a triangulation to large ! "
                "For machine's safety Simulator is going to terminate !");
    
    if (pbc)
        throw NotImplementedException("NotImplementedException:\n"
                "Periodic boundary conditions are  not implemented yet."
                "Simulator is about to terminate !");
}

void Simulator::addCell(const Cell& newCell)
{
    try
    {
        if (cells.size() == MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached."
                                    "New cell will not be added !");
        
        cells.push_back(newCell);
        numberofCells++;        
    } 
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
}

void Simulator::addCell()
{
    
    try
    {
        if (cells.size() == MAX_CELLS) 
            throw MaxSizeException("Maximum number of cells reached."
                                    "New cell will not be added !");
        
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
        cout << e.what() << endl;
        return;
    }
}

void Simulator::calcForces()
{
    // RESET FORCES
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].voidForces();
    }
    
    // CALCULATE INTRA-CELLULAR FORCES 
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces();
    }
    
    // CALCULATE FORCES BETWEEN CELLS AND BOX
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces(bs);
    }    
    
    // CALCULATE INTER-CELLULAR FORCES
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

void Simulator::doStep()
{
    calcForces();
    integrate();
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
    renderScript(drawBox);
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

void Simulator::renderScript(bool box)
{
    ofstream os;
    os.open(script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    if (box)
    {
        printBox(os);
    }
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
    os << "cmd.do(\"bg white\")\n\n";
    
    if (box)
    {
        os << "B = Box( ("<< -bs << "," << bs <<"), ("<< -bs << "," << bs <<"), ("<< -bs << "," << bs <<"), ";
        os << 2*bs << ", 2.5, color=(0.0,0.0,0.0) )\n";
        os << "obj = B.box\n";
        os << "cmd.load_cgo(obj, \"box\", 1)\n";
    }
    
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


void Simulator::printBox(ofstream& os)
{
    os << "class Box(object):\n";
    os << "  def __init__ (self, x_ax, y_ax, z_ax, step, linewidth, color):\n";
    os << "    lw = linewidth\n";
    os << "    c1 = color[0]\n";
    os << "    c2 = color[1]\n";
    os << "    c3 = color[2]\n";
    os << "\n";
    os << "    self.box=[]  \n";   
    os << "\n";
    os << "    self.box.append(LINEWIDTH)\n";
    os << "    self.box.append(float(lw))\n";
    os << "    self.box.append(BEGIN)\n";
    os << "    self.box.append(LINES)\n";
    os << "    self.box.append(COLOR)\n";
    os << "    self.box.append(c1)\n";
    os << "    self.box.append(c2)\n";
    os << "    self.box.append(c3)\n";
    os << "\n";  
    os << "    X = range(x_ax[0],x_ax[1]+1,step)\n";
    os << "    Y = range(y_ax[0],y_ax[1]+1,step)\n";
    os << "    Z = range(z_ax[0],z_ax[1]+1,step)\n"; 
    os << "\n";
    os << "    x_end = x_ax[1] - x_ax[0]\n";
    os << "    y_end = y_ax[1] - y_ax[0]\n";
    os << "    z_end = z_ax[1] - z_ax[0]\n";
    os << "\n";
    os << "    l = len(X)\n"; 
    os << "\n";    
    os << "    for i in range(len(X)):\n";
    os << "        for j in range(len(Y)):\n";
    os << "            x1 = X[i]\n";
    os << "            y1 = Y[j]\n";
    os << "            z1 = z_ax[0]\n";
    os << "            x2 = x1\n";
    os << "            y2 = y1\n";
    os << "            z2 = z1 + z_end\n";
    os << "\n";    
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x1))\n";
    os << "            self.box.append(float(y1))\n";
    os << "            self.box.append(float(z1))\n";
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x2))\n";
    os << "            self.box.append(float(y2))\n";
    os << "            self.box.append(float(z2))\n";
    os << "\n";
    os << "    for i in range(len(X)):\n";
    os << "        for j in range(len(Z)):\n";
    os << "            x1 = X[i]\n";
    os << "            y1 = y_ax[0]\n";
    os << "            z1 = Z[j]\n";
    os << "            x2 = x1\n";
    os << "            y2 = y1 + y_end\n";
    os << "            z2 = z1\n";
    os << "\n";
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x1))\n";
    os << "            self.box.append(float(y1))\n";
    os << "            self.box.append(float(z1))\n";
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x2))\n";
    os << "            self.box.append(float(y2))\n";
    os << "            self.box.append(float(z2))\n";
    os << "\n";
    os << "    for i in range(len(Y)):\n";
    os << "        for j in range(len(Z)):\n";
    os << "            x1 = x_ax[0]\n";
    os << "            y1 = Y[i]\n";
    os << "            z1 = Z[j]\n";
    os << "            x2 = x1 + x_end\n";
    os << "            y2 = y1\n";
    os << "            z2 = z1\n";
    os << "\n";
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x1))\n";
    os << "            self.box.append(float(y1))\n";
    os << "            self.box.append(float(z1))\n";
    os << "            self.box.append(VERTEX)\n";
    os << "            self.box.append(float(x2))\n";
    os << "            self.box.append(float(y2))\n";
    os << "            self.box.append(float(z2))\n";
    os << "\n";
    os << "    self.box.append(END)\n\n\n";
}