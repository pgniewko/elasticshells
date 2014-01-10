#ifndef SIMULATOR_H
#define	SIMULATOR_H

class Simulator {
public:
    Simulator();
    Simulator(const Simulator& orig);
    virtual ~Simulator();
    
    int generate_pos();
    int read_pos(const char* filename);
    void write_pos(const char* filename, bool wrap=1);
private:

};

#endif	/* SIMULATOR_H */

