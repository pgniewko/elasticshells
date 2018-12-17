#include "Simulator.h"

utils::Logger Simulator::simulator_logs("simulator");
unsigned long Simulator::FORCE_EVALUATION_COUTER(0);

bool Simulator::RESTART_FLAG(false);

Simulator::Simulator(const arguments& args) : number_of_shells(0), box(0, 0, 0),
    traj(args.traj_file, args.box_file), log_sim(args.output_file, args.ob_config_file),
    restarter(args.topology_file, args.lf_file)
{
    try
    {
        diagnoseParams(args);
    }
    catch (NotImplementedException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    catch (NotAllowedException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    catch (DataException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    //params.log_step = args.log_step;
    params.d = args.d;
    //params.nbhandler = args.nb_flag;
    params.E_shell = args.E_shell;
    params.nu = args.nu;
    params.th = args.thickness;
    params.dt = args.dt;
    params.dp = args.dp;
    params.ddp = args.ddp;
    params.ttime = args.ttime;
    params.r_vertex = args.r_vertex;
    //params.draw_box = args.draw_box;
    params.const_volume = args.const_volume;
    params.nsteps = args.nsteps ? args.nsteps : (int)params.ttime / params.dt;
    params.platotype = args.platotype;

    integrator = new Integrator(this);

    setTriangulator(args.tritype);
    box.set_x(args.bsx);
    box.set_y(args.bsy);
    box.set_z(args.bsz);
    box.set_x_min(args.bsxe);
    box.set_y_min(args.bsye);
    box.set_z_min(args.bsze);
    box.set_pbc(args.pbc);
    box.set_E(args.E_wall);
    box.set_nu(args.nu);
    box.set_default_schedule(params.nsteps, 1, args.bsdx, args.bsdy, args.bsdz, 0.0, 0.0, 0.0);
    box.configure_scheduler(args.sch_config_file);
    OsmoticForce::setVolumeFlag(args.osmotic_flag);
    OsmoticForce::setEpsilon(args.eps);
    Shell::bending = args.bending;
    logParams();

    fc = ForcesCalculator(estimate_m(), args.pbc, args.bending);
}

Simulator::~Simulator()
{
    delete integrator;
}

void Simulator::diagnoseParams(arguments args)
{
    if (args.d < 0)
        throw NotImplementedException("NotImplementedException:\n"
                                      "Single point representation is not implemented yet. "
                                      "Simulator is about to terminate !");

    if ( (args.d == 0 && STRCMP(args.tritype, "rnd")) || (args.d == 0 && STRCMP(args.tritype, "plato")) )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Random and Plato triangulation cannot be used "
                                  "with single point representation.\n"
                                  "Simulation exits with EXIT_FAILURE status !");

    if ( (args.d == 0 && args.restart) || (args.d == 0 && args.analyze) )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Simulation for a single node representation"
                                  "cannot be restarted or post-analyzed!");

    if ( args.d == 0 && args.const_volume && args.nu != 0.5 )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Single node representation:"
                                  " for const-volume option nu has to be equal to 0.5!");

    if (args.d > 8)
        throw DataException("DataException:\n"
                            "Depth of a triangulation too large ! "
                            "For machine's safety Simulator is about to terminate !");

    if (args.d < 0)
        throw DataException("Depth of a triangulation cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.E_wall < 0)
        throw DataException("Effective shell-box Young's modulus cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.E_shell < 0)
        throw DataException("Effective shell-shell Young's modulus cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.dt <= 0)
        throw DataException("Time step must be positive number ! \n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.r_vertex <= 0)
        throw DataException("Vertex radius must be larger than 0! \n!"
                            "Simulation will terminate with exit(1)!\n");

}

void Simulator::logParams()
{
    simulator_logs << utils::LogLevel::INFO  << "SIM_STEPS=" << params.nsteps << "\n";
//    simulator_logs << utils::LogLevel::INFO  << "LOG_STEP="  << params.log_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "TRIANGULATOR="  << triangulator << "\n";
    simulator_logs << utils::LogLevel::FINE  << "TIME STEP(DT)="  << params.dt << " [s]\n";
    simulator_logs << utils::LogLevel::FINE  << "DEPTH="  << params.d << "\n";
    simulator_logs << utils::LogLevel::FINE  << "DP="  << params.dp << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "DDP="  << params.ddp << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "E SHELL="  << params.E_shell << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "E BOX="  << box.get_E() << " [MPa]\n";
    //simulator_logs << utils::LogLevel::FINE  << "SURFACE_MODULUS="  << (params.E_shell * params.th) << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (SHELL)="  << params.nu << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (BOX)="  << box.get_nu() << "\n";
    simulator_logs << utils::LogLevel::FINE  << "R_VERTEX="  << params.r_vertex << " [micron]\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.PBC=" << (box.pbc ? "true" : "false") << "\n";
    //simulator_logs << utils::LogLevel::FINE  << "BOX.BOX_DRAW=" << (params.draw_box ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_FLAG=" << (OsmoticForce::getFlag() ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_EPS=" << OsmoticForce::getEpsilon() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.X="  << box.get_x() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Y="  << box.get_y() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Z="  << box.get_z() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.XE=" << box.get_x_min() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.YE=" << box.get_y_min() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.ZE=" << box.get_z_min() << "\n";
}

void Simulator::init_shells(int N, double r_min, double r_max, bool jam)
{
    if (r_min > r_max)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Illegal arguments: r_min > r_max. Simulator sets: r_min = r_max \n";
        r_max = r_min;
    }

    simulator_logs << utils::LogLevel::INFO  << "BENDING: " << (!Shell::bending ? "true" : "false") << "\n";

    double nx, ny, nz;
    bool flag = true;
    double rc = 2.0 * params.r_vertex;
    double rsum;
    double r0;

    while (number_of_shells < N)
    {
        r0 = uniform(r_min, r_max);
        flag = true;
        nx = uniform(-box.get_x() + r0 + rc, box.get_x() - r0 - rc);
        ny = uniform(-box.get_y() + r0 + rc, box.get_y() - r0 - rc);
        nz = uniform(-box.get_z() + r0 + rc, box.get_z() - r0 - rc);

        if (number_of_shells == 0)
        {
            nx = 0;
            ny = 0;
            nz = 0;
        }

        Vector3D shift(nx, ny, nz);

        for (int i = 0; i < number_of_shells; i++)
        {
            Vector3D tmpcm = shells[i].get_cm();
            Vector3D delta = shift - tmpcm;
            rsum = shells[i].get_r0() + r0 + rc + constants::epsilon;

            if ( delta.length() < rsum)
            {
                flag = false;
            }
        }

        if (flag)
        {
            addShell(r0);
            shiftShell(shift, number_of_shells - 1);
        }
    }

    if (jam)
    {
        simulator_logs << utils::LogLevel::INFO  << "SIMULATION STARTS FROM THE JAMMED PACKING\n";

        if (params.d == 0)
        {
            Packer::packShells(box, shells, params.th, false);
        }
        else
        {
            Packer::packShells(box, shells, params.th, true);
        }

        fc.reset_dl(estimate_m(), box);
    }

    if (params.d == 0)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Single vertex representation\n";
        simulator_logs << utils::LogLevel::WARNING  << "Vertex radius reassignment: shell.vertex_r = shell.init_r\n";

        for (int i = 0; i < number_of_shells; i++)
        {
            shells[i].set_vertex_size( shells[i].get_r0() );
        }
    }

    restarter.saveTopologyFile(shells);
    set_min_force();

    create_shells_image();
    copy_shells_data();
}

void Simulator::pushShell(const Shell& newShell)
{
    try
    {
        shells.push_back(newShell);
        number_of_shells = shells.size();
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::WARNING << e.what() << "\n";
        return;
    }
}

void Simulator::addShell(double r0)
{
    try
    {
        std::list<Triangle> tris;

        if ( !triangulator.compare("plato") )
        {
            PlatonicTriangulatoin tio(params.d, params.platotype);
            tris = tio.triangulate(r0);
        }
        else if ( !triangulator.compare("simple") )
        {
            SimpleTriangulation sm(params.d);
            tris = sm.triangulate(r0);
        }

        else if ( !triangulator.compare("rnd") )
        {
            RandomTriangulation rnd(10, 100, 0.1, 1000.0, params.r_vertex);
            tris = rnd.triangulate(r0);
        }

        Shell new_shell(tris);
        new_shell.set_vertex_size(params.r_vertex);
        new_shell.set_shell_id(number_of_shells);
        new_shell.set_r0(r0);

        new_shell.set_ecc(params.E_shell);
        new_shell.set_nu(params.nu);
        new_shell.set_elements_parameters(params.E_shell, params.th, params.nu);
        new_shell.set_hinges(params.E_shell, params.th, params.nu);
        new_shell.set_dp(params.dp, params.ddp);

        double radial_eps = 1.0 + (0.5 * (1 - params.nu)  * (new_shell.get_turgor() * r0) / (params.E_shell * params.th));

        if (new_shell.get_number_vertices() == 1)
        {
            radial_eps = 1.0;
        }

        new_shell.set_constant_volume (radial_eps);

        pushShell(new_shell);
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}

void Simulator::restart()
{
    simulator_logs << utils::LogLevel::INFO << "Simulation runs in [restart] mode. \n" ;
    restarter.registerVMap();
    restarter.readTopologyFile(shells);
    restarter.readLastFrame(shells);
    number_of_shells = shells.size();
    set_min_force();
    Simulator::RESTART_FLAG = true;
    create_shells_image();
    copy_shells_data();
}

void Simulator::analyze()
{
    restarter.registerVMap();
    restarter.readTopologyFile(shells);
    number_of_shells = shells.size();

    uint frames_number = traj.count_frames();
    simulator_logs << utils::LogLevel::INFO << " Number of frames in a trajectory file: " << (int) frames_number << "\n" ;

    std::vector<std::string> turgor_list = log_sim.read_turgors_file();
    std::vector<std::string> boxsize_list = traj.read_saved_box();

    if (turgor_list.size() != frames_number)
    {
        simulator_logs << utils::LogLevel::SEVERE  << "CORRUPTED DATA: Turgor data doesn't match the number of frames. Observation analysis exists!";
        exit(EXIT_FAILURE);
    }

    if (boxsize_list.size() != frames_number)
    {
        simulator_logs << utils::LogLevel::SEVERE  << "CORRUPTED DATA: Box size data doesn't match the number of frames. Observation analysis exists!";
        exit(EXIT_FAILURE);
    }

    log_sim.register_observers();
    log_sim.open();
    log_sim.print_header();

    for (std::size_t i = 1; i <= frames_number; i++)
    {
        simulator_logs << utils::LogLevel::INFO << "[analyze] Processing frame number: " << i <<  "/" <<  (int) frames_number << "\n" ;
        restarter.read_frame(traj.get_traj_file(), shells, i);
        restarter.assign_turgors(turgor_list[i - 1], shells);
        restarter.assign_box_size(boxsize_list[i - 1], box);
        recalculate_mass_centers();

        if (i == 1)
        {
            set_min_force();
            simulator_logs << utils::LogLevel::FINE  << "MIN_FORCE (in <<analyze>> mode) SET TO= "  << sqrt(MIN_FORCE) << " [units?]\n";
        }

        create_shells_image();
        copy_shells_data();
        log_sim.dump_state(box, shells);
    }
}

void Simulator::simulate(int steps)
{
    // LOGGER READY TO WORK
    log_sim.register_observers();
    log_sim.open();
    log_sim.print_header();

    // TRAJECTORY FILE OPEND FOR DUMP
    traj.open();

    // IF SIMULATION RESTART - DON'T DUMP THE STATS
    if (!Simulator::RESTART_FLAG)
    {
        traj.save_traj(shells, getTotalVertices());
        log_sim.dump_state(box, shells);
        saveTurgors();
        restarter.saveLastFrame(shells, box);
        restarter.saveTopologyFile(shells);
        traj.save_box(box, steps * params.dt);
        box.save_remaining_schedule();
    }

    calcForces();

    bool resized = false;

    for (int step = 0; step < steps; step++)
    {
        if ( box.nthTodo() )
        {
            simulator_logs << utils::LogLevel::INFO << "The simulation has reached the end of the schedule @ the step " << step << "/" << steps << ".\n";
            break;
        }

        if ( step % (steps / std::min(steps, 10) ) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * step / steps << "% OF THE SIMULATION IS DONE" "\n";
        }

        copy_shells_data();

        do
        {
            do
            {
                integrate();
            }
            while ( check_min_force() );
        }
        while ( check_const_volume() );

        copy_back_shells_data();
        integrator->resetParams(this);


        //if ( (step + 1) % params.log_step == 0 )
        //{
            // ** SAVE COORDINATES - i.e. "logging" coordinates
            traj.save_traj(shells, getTotalVertices());
            log_sim.dump_state(box, shells);
            saveTurgors();
            traj.save_box(box, (step + 1) * params.dt);
            restarter.saveLastFrame(shells, box);
            restarter.saveTopologyFile(shells);
        //}

        if ( step < steps - 1 ) // DO NOT RESIZE ON THE LAST STEP
        {
            resized = box.resize( volumeFraction() );
            fc.reset_dl(estimate_m(), box);
        }
        else
        {
            resized = false;
        }

        if (resized)
        {
            box.save_remaining_schedule();
        }

        recenterShells();
    }

    log_sim.dump_state(box, shells);
    saveTurgors();
    traj.save_box(box, steps * params.dt);
    restarter.saveLastFrame(shells, box);
    traj.save_traj(shells, getTotalVertices());
    box.save_remaining_schedule();

    traj.close();
    log_sim.close();

    simulator_logs << utils::LogLevel::FINEST << "Forces have been evaluated " << FORCE_EVALUATION_COUTER << " times.\n";
    simulator_logs << utils::LogLevel::FINEST << "Energy has been evaluated "  << Energy::ENERGY_EVALUATION_COUNTER << " times.\n";
}

void Simulator::calcForces()
{
    for (uint i = 0; i < forces.size(); i++)
    {
        forces[i] = 0.0;
    }

    fc.calculate_forces(xyz, forces, elements, hinges, vs_map, graph, turgors, shells.size(),
                        params.r_vertex, params.E_shell, params.nu,
                        box.get_E(), box.get_nu());

    FORCE_EVALUATION_COUTER++;

}

void Simulator::shiftShell(const Vector3D& v3d, int shell_id)
{
    shells[shell_id].add_vector(v3d);
    shells[shell_id].calc_cm(); //.update();
}

double Simulator::volumeFraction()
{
    double box_vol = box.get_volume();
    double shells_volume = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        shells_volume += shells[i].calc_volume();
    }

    return (shells_volume / box_vol);

}

int Simulator::getTotalVertices()
{
    int totalnumber = 0;

    for (int i = 0; i < number_of_shells; i++)
    {
        totalnumber += shells[i].get_number_vertices();
    }

    return totalnumber;
}

double Simulator::getLengthScale(double r_0)
{
    double maxscale = 0.0;

    maxscale = std::max(maxscale, params.r_vertex);

    if (params.d == 0)
    {
        maxscale = std::max(maxscale, r_0) ;
    }

    return maxscale;
}

/*
 * TRIANGULATION
 */
void Simulator::setTriangulator(char* token)
{
    if (STRCMP (token, "simple"))
    {
        triangulator = std::string(token);
    }
    else if (STRCMP (token, "plato"))
    {
        triangulator = std::string(token);
    }
    else if (STRCMP (token, "rnd"))
    {
        triangulator = std::string(token);
    }
    else
    {
        triangulator = std::string("simple");
        simulator_logs << utils::LogLevel::FINE  << "SIMPLE TRIANGULATION IS APPLIED\n";
    }
}

void Simulator::recalculate_mass_centers()
{
    for (uint i = 0; i < shells.size(); i++)
    {
        shells[i].calc_cm();
    }
}

void Simulator::set_min_force()
{
    double average_area = shells[0].calc_surface_area();
    average_area /= shells[0].get_number_triangles();

    double max_turgor = 0.0;

    for (int i = 0; i < number_of_shells; i++)
    {
        max_turgor = std::max(max_turgor, shells[i].get_turgor());
    }

    MIN_FORCE = FORCE_FRAC * max_turgor * average_area;

    if (shells[0].get_number_vertices() == 1)
    {
        MIN_FORCE = 1e-12;
    }

    Shell::FORCE_FRAC   = FORCE_FRAC;
    Shell::MIN_FORCE    = MIN_FORCE;

    simulator_logs << utils::LogLevel::FINE  << "MIN_FORCE = "  << sqrt(MIN_FORCE) << " [units?]\n";
}

bool Simulator::check_min_force()
{
    for (uint i = 0; i < forces.size(); i++)
    {
        if (constants::sqrt3 * forces[i] > MIN_FORCE)
        {
            return true;
        }
    }

    return false;
}

bool Simulator::check_const_volume()
{
    if (params.const_volume)
    {
        copy_back_shells_data();
        double step;
        double eps = 0.001; // 0.1% accuracy
        bool flag = false;

        for (int i = 0; i < number_of_shells; i++)
        {
            step = shells[i].check_volume_condition();

            if ( fabs(step) > eps )
            {
                flag = true;
                turgors[i] = shells[i].ajust_turgor(step);
            }
        }

        return flag;
    }
    else
    {
        return false;
    }
}

void Simulator::integrate()
{
    integrator->integrate(this);
}

void Simulator::saveTurgors()
{
    std::string turgorDumpFile = log_sim.get_file_name() + ".turgor.out";

    double box_volume = box.get_volume();
    double box_x = box.get_x();
    double box_y = box.get_y();
    double box_z = box.get_z();
    double shells_volume = 0.0;

    FILE* os = fopen(turgorDumpFile.c_str(), "a");

    if ( os == NULL )
    {
        os = fopen(turgorDumpFile.c_str(), "w");
    }

    if ( os == NULL )
    {
        simulator_logs << utils::LogLevel::WARNING << "Can not open file:<<" << turgorDumpFile << "for writing.\n";
    }

    double turgor;
    Vector3D cm;

    for (uint i = 0; i < shells.size(); i++)
    {
        shells_volume += shells[i].calc_volume();
    }

    fprintf(os, "%-7.4f %5.3f %5.3f %5.3f %7.4f ", box_volume, box_x, box_y, box_z, shells_volume);


    for (uint i = 0; i < shells.size(); i++)
    {
        turgor = shells[i].get_turgor();
        cm = shells[i].get_cm();

        fprintf(os, "%5.4f %6.4f %6.4f %6.4f ", turgor, cm.x, cm.y, cm.z);
    }

    fprintf(os, "%s", "\n");

    if ( os != NULL )
    {
        fclose(os);
    }
}

void Simulator::recenterShells()
{
    if (box.pbc)
    {
        for (uint i = 0; i < shells.size(); i++)
        {
            Vector3D shift = Box::recentering_vector( shells[i].get_cm(), box );
            shells[i].add_vector(shift);
        }
    }

    return;
}

void Simulator::create_shells_image()
{
    simulator_logs << utils::LogLevel::INFO  << "Creating local image of the shells\n";
    int vertex_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        turgors.push_back( 0.0 );
        double x_, y_, z_;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            x_ = shells[i].vertices[j].r_c.x;
            y_ = shells[i].vertices[j].r_c.y;
            z_ = shells[i].vertices[j].r_c.z;

            xyz.push_back(x_);
            xyz.push_back(y_);
            xyz.push_back(z_);

            forces.push_back(0.0);
            forces.push_back(0.0);
            forces.push_back(0.0);

            object_map vm(i, j);
            vs_map.push_back(vm);

            inv_vs_map[vm] = vertex_counter;

            std::vector<int> bonds;
            graph.push_back(bonds);

            vertex_counter++;
        }
    }

    for (uint i = 0; i < shells.size(); i++)
    {
        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            object_map vm_ij(i, j);
            int ij_id = inv_vs_map[vm_ij];

            for (int k = 0; k < shells[i].get_number_vertices(); k++)
            {
                if (k != j && shells[i].vertices[j].is_neighbor(k))
                {
                    object_map vm_ik(i, k);
                    int ik_id = inv_vs_map[vm_ik];
                    graph[ij_id].push_back(ik_id);
                }
            }

        }
    }


    int element_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        int ia, ib, ic;
        int ia_mapped, ib_mapped, ic_mapped;

        for (int j = 0; j < shells[i].get_number_triangles(); j++)
        {
            element el_;
            ia = shells[i].triangles[j].ia;
            ib = shells[i].triangles[j].ib;
            ic = shells[i].triangles[j].ic;

            object_map vm_a(i, ia);
            object_map vm_b(i, ib);
            object_map vm_c(i, ic);

            ia_mapped = inv_vs_map[vm_a];
            ib_mapped = inv_vs_map[vm_b];
            ic_mapped = inv_vs_map[vm_c];

            el_.ia = ia_mapped;
            el_.ib = ib_mapped;
            el_.ic = ic_mapped;

            el_.an[0] = shells[i].triangles[j].an[0];
            el_.an[1] = shells[i].triangles[j].an[1];
            el_.an[2] = shells[i].triangles[j].an[2];

            el_.L2[0] = shells[i].triangles[j].L2[0];
            el_.L2[1] = shells[i].triangles[j].L2[1];
            el_.L2[2] = shells[i].triangles[j].L2[2];

            el_.ki[0] = shells[i].triangles[j].ki[0];
            el_.ki[1] = shells[i].triangles[j].ki[1];
            el_.ki[2] = shells[i].triangles[j].ki[2];

            el_.ci[0] = shells[i].triangles[j].ci[0];
            el_.ci[1] = shells[i].triangles[j].ci[1];
            el_.ci[2] = shells[i].triangles[j].ci[2];

            elements.push_back(el_);

            object_map vm(i, j);
            ts_map.push_back(vm);
            inv_ts_map[vm] = element_counter;

            element_counter++;
        }
    }

    int hinge_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        int x1, x2, x3, x4;
        int x1_mapped, x2_mapped, x3_mapped, x4_mapped;

        for (int j = 0; j < shells[i].get_number_hinges(); j++)
        {

            hinge h_;
            x1 = shells[i].hinges[j].x1;
            x2 = shells[i].hinges[j].x2;
            x3 = shells[i].hinges[j].x3;
            x4 = shells[i].hinges[j].x4;

            object_map vm_a(i, x1);
            object_map vm_b(i, x2);
            object_map vm_c(i, x3);
            object_map vm_d(i, x4);

            x1_mapped = inv_vs_map[vm_a];
            x2_mapped = inv_vs_map[vm_b];
            x3_mapped = inv_vs_map[vm_c];
            x4_mapped = inv_vs_map[vm_d];

            h_.v1 = x1_mapped;
            h_.v2 = x2_mapped;
            h_.v3 = x3_mapped;
            h_.v4 = x4_mapped;

            h_.D = shells[i].hinges[j].D;
            h_.theta0 = shells[i].hinges[j].theta0;
            h_.sinTheta0 = shells[i].hinges[j].sinTheta0;

            hinges.push_back(h_);


            object_map vm(i, j);
            hs_map.push_back(vm);
            inv_hs_map[vm] = hinge_counter;
            hinge_counter++;
        }
    }

    integrator->set_n( xyz.size() );
}

void Simulator::copy_shells_data()
{
    int vertex_no = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        double x_, y_, z_;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            x_ = shells[i].vertices[j].r_c.x;
            y_ = shells[i].vertices[j].r_c.y;
            z_ = shells[i].vertices[j].r_c.z;

            object_map vm(i, j);

            vertex_no =  inv_vs_map[ vm ];

            xyz[3 * vertex_no + 0] = x_;
            xyz[3 * vertex_no + 1] = y_;
            xyz[3 * vertex_no + 2] = z_;
        }
    }

    for (uint i = 0; i < shells.size(); i++)
    {
        turgors[i] = shells[i].get_turgor();
    }

    fc.set_dl_dims(-box.get_x(), box.get_x(), 0);
    fc.set_dl_dims(-box.get_y(), box.get_y(), 1);
    fc.set_dl_dims(-box.get_z(), box.get_z(), 2);
}


void Simulator::copy_back_shells_data()
{
    for (uint i = 0; i < vs_map.size(); i++)
    {
        int shell_id = vs_map[i].shell_id;
        int vert_id  = vs_map[i].vert_id;
        shells[shell_id].vertices[vert_id].r_c.x = xyz[3 * i + 0];
        shells[shell_id].vertices[vert_id].r_c.y = xyz[3 * i + 1];
        shells[shell_id].vertices[vert_id].r_c.z = xyz[3 * i + 2];
    }
}

int Simulator::estimate_m()
{
    double x_dim = 2 * box.get_x();
    double y_dim = 2 * box.get_y();
    double z_dim = 2 * box.get_z();

    double d = 2 * params.r_vertex;

    int m = std::min( (int)floor(x_dim / d), (int) floor(y_dim / d) );
    m = std::min(m, (int) floor(z_dim / d) );
    m = std::min(m, MAX_M);
    return m;
}