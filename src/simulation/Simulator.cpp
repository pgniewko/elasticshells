#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args), numberofCells(0), box(0,0,0), 
        logStep(params.log_step), saveStep(params.save_step), boxStep(params.box_step), vlistStep(params.vlist_step)
{
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
//    bs = params.bs;
    drawBox = params.draw_box;
    pbc = params.pbc;
    gamma = params.k;
    Rc = params.r_cut;
    verlet_r = params.verlet_r;
    ttotal = params.ttime;
    trajfile = params.traj_file;
    script = params.output_file;
    surfaceScript = params.surface_file;
    nsteps = (int) ttotal / dt;
    setIntegrator(params.integrator_a);
    
    box.setX(params.bsx);
    box.setY(params.bsy);
    box.setZ(params.bsz);
    box.setDx(params.bsdx);
    box.setDy(params.bsdy);
    box.setDz(params.bsdz);
    
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

Simulator::Simulator(const Simulator& orig) : box(orig.box){}

Simulator::~Simulator() {}

void Simulator::diagnoseParams()
{
    if (d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                "Single point representation is not implemented yet. "
                "Simulator is about to terminate !");
    if (d > 9)
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
    addCell(1.0);
}

void Simulator::addCell(double r0)
{
    
    try
    {
        if (cells.size() == MAX_CELLS) 
            throw MaxSizeException("Maximum number of cells reached."
                                    "New cell will not be added !");
        
        SimpleTriangulation sm(params.d);
        list<Triangle> tris = sm.triangulate(r0);
        Cell newCell(tris);
   
        newCell.setA(a);
        newCell.setDp(dp);
        newCell.setRc(Rc);
        newCell.setGamma(gamma);
        newCell.setVerletR(verlet_r);
        newCell.setCellId(numberofCells);
        newCell.setInitR(r0);
        
        addCell(newCell);       
    } 
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
}

void Simulator::initCells(int N, double r0)
{
    
    double nx, ny, nz;
    bool flag = true;
    while(numberofCells < N)
    {
        flag = true; 
        nx = uniform(-box.getX()+r0+Rc, box.getX()-r0-Rc);
        ny = uniform(-box.getY()+r0+Rc, box.getY()-r0-Rc);
        nz = uniform(-box.getZ()+r0+Rc, box.getZ()-r0-Rc);
        Vector3D shift(nx, ny, nz);
        for (int i = 0; i < numberofCells; i++)
        {
            Vector3D tmpcm = cells[i].getCm();
            Vector3D delta = shift - tmpcm;
            double rsum = cells[i].getInitR() + r0; 
            if ( (delta.length() + EPSILON) < rsum)
            {
               flag = false; 
            }
        }
        if (flag)
        {
            addCell(r0);
            moveCell(shift, numberofCells-1);
        }

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
        cells[i].calcForces(box);
    }    
    
    // CALCULATE INTER-CELLULAR FORCES
    for (int i = 0; i < numberofCells; i++) 
    {
        for (int j = 0; j < numberofCells; j++) 
        {
            //if (i != j)
            //{
                //cells[i].calcForces(cells[j]);
                cells[i].calcForcesVL(cells[j]);   
            //}
        }
    }
}

void Simulator::integrateEuler()
{
    calcForces();
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / m;
        }
    
    }
}

void Simulator::heunMethod()
{
    calcForces();
    double m;
    
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
            cells[i].vertices[j].tmp_force = cells[i].vertices[j].force;
        }
    }
    
    //move the whole time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / m;
        }
    }
    calcForces();
    
    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + 0.5 * dt * ( cells[i].vertices[j].tmp_force + cells[i].vertices[j].force) / m;
        }
    }  
}

void Simulator::midpointRungeKutta()
{
    calcForces();
    double m;
    
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
        }
    }
    
    //move half time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += 0.5 * dt * cells[i].vertices[j].force / m;
        }
    }
    calcForces();
    
    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + dt * cells[i].vertices[j].force / m;
        }
    }    
}

void Simulator::integrateVv()
{
    calcForces();
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
        }
    }
    
    calcForces();
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; // v(t+1) = v(t+1)_1 + 0.5 * dt * a(t+1)
        }
    }
}

void Simulator::simulate()
{
    simulate(nsteps);       
}

void Simulator::simulate(int steps)
{   
    rebuildVerletLists();
    renderScript(drawBox);
    saveSurfaceScript();
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
        lastCellIndex += cells[i].numberofVertices();
    }

    for (int i = 0; i < steps; i++)
    {   
        if ( (i+1) % boxStep == 0)
        {
            box.resize();
            cout << "i= " << i <<" " << box.getX() << " " << box.getY() << " " << box.getZ() << " " << endl;
        }
        
        
        integrate();
        if ( (i+1) % vlistStep == 0)
        {
            rebuildVerletLists();
        }
        
        if ( (i+1) % saveStep == 0)
        {
            
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
                lastCellIndex += cells[i].numberofVertices();
            }
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
    if (STRCMP (token, "fe"))
    {
        this->setIntegrator(&Simulator::integrateEuler);
    }
    else if (STRCMP (token, "hm"))
    {
        this->setIntegrator(&Simulator::heunMethod);
    }
    else if (STRCMP (token, "rk"))
    {
        this->setIntegrator(&Simulator::midpointRungeKutta);
    }
    else if (STRCMP (token, "vv"))
    {
        this->setIntegrator(&Simulator::integrateVv);
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

void Simulator::rebuildVerletLists()
{
    for (int i = 0; i < numberofCells; i++)
    {
        cells[i].voidVerletLsit();
    }
    
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < numberofCells; j ++)
        {
           cells[i].builtVerletList(cells[j]);    
        }
    }
}

void Simulator::renderScript(bool boxFlag)
{
    ofstream os;
    os.open(script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    if (boxFlag)
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
        lastCellIndex += cells[i].numberofVertices();
    }
    
    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n\n";
    
    if (boxFlag)
    {
        os << "B = Box(";
        os << "("<< -box.getX() << "," << box.getX() <<"),";
        os << "("<< -box.getY() << "," << box.getY() <<"),";
        os << "("<< -box.getZ() << "," << box.getZ() <<"),";
        os << " 2.5, color=(0.0, 0.0, 0.0) )\n";
        os << "obj = B.box\n";
        os << "cmd.load_cgo(obj, \"box\", 1)\n";
    }
    
    os.close();
}

void Simulator::saveSurfaceScript()
{
    ofstream os;
    os.open(surfaceScript);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    

    os << "def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):\n\n";

    os << "  nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)\n";
    os << "  ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)\n";
    os << "  nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)\n\n";

    os << "  return (nx,ny,nz)\n\n\n";

    os << "def draw_plane_cgo(name, apex1, apex2, apex3, color=(1,1,1)):\n";
    os << "  x1,y1,z1 = map(float,apex1)\n";
    os << "  x2,y2,z2 = map(float,apex2)\n";
    os << "  x3,y3,z3 = map(float,apex3)\n";
    os << "  if type(color) == type(''):\n";
    os << "    color = map(float,color.replace('(','').replace(')','').split(','))\n\n";

    os << "  # Compute the normal vector for the triangle\n";
    os << "  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)\n\n";

    os << "  # Create the CGO objects\n";
    os << "  obj = [\n";

    os << "    BEGIN, TRIANGLE_STRIP,\n\n";

    os << "    COLOR, color[0], color[1], color[2],\n";
    os << "    NORMAL, normal1[0], normal1[1], normal1[2],\n";
    os << "    NORMAL, -normal1[0], -normal1[1], -normal1[2],\n";
    os << "    VERTEX, x1, y1, z1,\n";
    os << "    VERTEX, x2, y2, z2,\n";
    os << "    VERTEX, x3, y3, z3,\n\n";

    os << "    END\n";
    os << "  ]\n\n";

    os << "  # Display them\n";
    os << "  cmd.load_cgo(obj, name)\n\n";

    os << "def draw_plane(name,atom1='(pk1)',atom2='(pk2)',atom3='(pk3)',color=(1,1,1)):\n";
    os << "# get coordinates from atom selections\n";
    os << "  coor1 = cmd.get_model(atom1).atom[0].coord\n";
    os << "  coor2 = cmd.get_model(atom2).atom[0].coord\n";
    os << "  coor3 = cmd.get_model(atom3).atom[0].coord\n";
    os << "  draw_plane_cgo(name,coor1,coor2,coor3,color)\n\n";

    os << "cmd.extend(\"draw_plane\", draw_plane)\n\n\n";
          
    int iidx, jidx;
    int lastCellIndex = 0;
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            iidx = cells[i].vertices[j].getId() + 1 + lastCellIndex;
            os << "cmd.do(\"select H "<< iidx << ", name H" << iidx << "\")\n";
        }
        lastCellIndex += cells[i].numberofVertices();
    }
    
    
    lastCellIndex = 0;
    int faceCounter = 0;
    int idxa, idxb, idxc;
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberofFaces(); j++)
        {
            idxa = cells[i].triangles[j].ia + 1 + lastCellIndex;
            idxb = cells[i].triangles[j].ib + 1 + lastCellIndex;
            idxc = cells[i].triangles[j].ic + 1 + lastCellIndex;
            os << "cmd.do(\"draw_plane \\\"face"<<faceCounter<<"\\\", H_"<< idxa << ", H_" << idxb << ", H_" << idxc<< ", (0.8, 0.8, 0.8) \")\n";
            faceCounter++;
        }
        lastCellIndex += cells[i].numberofVertices();
    }
    
}

void Simulator::saveCellsState(const char* fileout)
{
    int index;
    ofstream os;
    os.open(fileout);
//    cout << fileout << endl;
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

//TODO: make it safe
Cell Simulator::getCell(int cell_index)
{
    return cells[cell_index];
}

void Simulator::saveCellsState()
{
    saveCellsState(params.output_file);
}

void Simulator::setBoxSize(const double bs)
{
    box.setX(bs);
    box.setY(bs);
    box.setZ(bs);
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
  os << "  def __init__ (self, x, y, z, linewidth, color):\n";
  os << "    lw = linewidth\n";
  os << "    c1 = color[0]\n";
  os << "    c2 = color[1]\n";
  os << "    c3 = color[2]\n";

  os << "    self.box = [\n";
  os << "    LINEWIDTH, float(lw), BEGIN, LINES,\n";
  os << "    COLOR,  color[0], color[1], color[2],\n";

  os << "    VERTEX, x[0], y[0], z[0],\n";
  os << "    VERTEX, x[1], y[0], z[0],\n";
  os << "    VERTEX, x[0], y[0], z[0],\n";
  os << "    VERTEX, x[0], y[1], z[0],\n";
  os << "    VERTEX, x[0], y[0], z[0],\n"; 
  os << "    VERTEX, x[0], y[0], z[1],\n"; 

  os << "    VERTEX, x[1], y[1], z[1],\n";
  os << "    VERTEX, x[1], y[1], z[0],\n"; 
  os << "    VERTEX, x[1], y[1], z[1],\n";
  os << "    VERTEX, x[0], y[1], z[1],\n";    
  os << "    VERTEX, x[1], y[1], z[1],\n";
  os << "    VERTEX, x[1], y[0], z[1],\n"; 

  os << "    VERTEX, x[0], y[0], z[1],\n"; 
  os << "    VERTEX, x[1], y[0], z[1],\n"; 

  os << "    VERTEX, x[0], y[1], z[0],\n"; 
  os << "    VERTEX, x[0], y[1], z[1],\n"; 

  os << "    VERTEX, x[0], y[1], z[0],\n"; 
  os << "    VERTEX, x[1], y[1], z[0],\n";

  os << "    VERTEX, x[0], y[0], z[1],\n"; 
  os << "    VERTEX, x[0], y[1], z[1],\n";      

  os << "    VERTEX, x[1], y[0], z[0],\n"; 
  os << "    VERTEX, x[1], y[1], z[0],\n";

  os << "    VERTEX, x[1], y[0], z[0],\n";
  os << "    VERTEX, x[1], y[0], z[1],\n";

  os << "    END\n";
  os << "    ]\n\n\n";
}