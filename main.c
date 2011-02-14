// MC Algorithm
// Michael Dacko, Lewin Stein, Christoph Tavan
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "ran.h"

#define PI 3.141592654

/**
 * Structure definitions
 */
struct TVector {
	double x;
	double y;

	TVector& operator=(TVector const& vector2) {
		if(this != &vector2) {
			x = vector2.x;
			y = vector2.y;
		}
		return *this;
	}

	TVector& operator=(double const& value) {
		x = value;
		y = value;
		return *this;
	}

	TVector& operator+=(TVector const& vector) {
		x += vector.x;
		y += vector.y;
		return *this;
	}

	TVector& operator-=(TVector const& vector) {
		x -= vector.x;
		y -= vector.y;
		return *this;
	}

	TVector& operator/=(double const& value) {
		x /= value;
		y /= value;
		return *this;
	}

	TVector& operator*=(double const& value) {
		x *= value;
		y *= value;
		return *this;
	}

	TVector operator*(double const& value) const {
		TVector temp(*this);
		temp.x *= value;
		temp.y *= value;
		return temp;
	}

	TVector operator/(double const& value) const {
		TVector temp(*this);
		temp.x /= value;
		temp.y /= value;
		return temp;
	}

	TVector operator+(TVector const& vector) const {
		TVector temp(*this);
		temp.x += vector.x;
		temp.y += vector.y;
		return temp;
	}

	TVector operator-(TVector const& vector) const {
		TVector temp(*this);
		temp.x -= vector.x;
		temp.y -= vector.y;
		return temp;
	}
};

struct TPartDist {				// for particle distances
	double dx;
	double dy;
	double r2;
};

struct TAverager
{
	double sum;
	double ssum;
	int number;
	
	void init()
	{
		sum = 0.0;
		ssum = 0.0;
		number = 0;
	}

	void add(double const& value)
	{
		sum += value;
		// ssum += value*value;
		number++;
	}

	double average()
	{
		return sum/(double)number;
	}
};

double rand_value(double min, double max);
TPartDist distance_special(double x1, double y1, double x2, double y2);
void print_coords(const char* filename);
void update_histogram();

void init_lj_shift();
double u_lj(double r2);
double u_lj_shifted(double r2);
double u_lj_deriv(double r2);
double u_lj_deriv_shifted(double r2);
TVector f_lj(TPartDist dr);
TVector f_lj_shifted(TPartDist dr);
double u_lj_shift;
double u_lj_deriv_shift;

void forces();
void scale_velocities();
void update_histogram();

/**
 * Global variables
 */
const bool debug = false;			// Enable/disable debug output
Ran myran(time(NULL));				// Random number generator (see ran.h)
const int nsamp = 10;				// Sample every nsamp timesteps
const double grmax = 3.5;			// Maximum radius for sampling g(r)

// System parameters
int N;						// Number of particles
double rho;					// Density of the system
double lx = 0, ly = 0;		// Dimensions of the simulation box, is calculated according 
							// to the density and Number of Particles.
double T;					// Temperature
double rc, rc2;				// LJ-cutoff-radius
TVector* r = NULL;			// Arrays for particle positions
TVector* r_next = NULL;
TVector* v = NULL;			// Arrays for particle velocities
TVector* v_next = NULL;
TVector* F = NULL;			// Arrays for inter-particle-forces

TVector vsum;				// Velocity of center of mass
double vsum2 = 0;			// Mean squared velocity = 2 E_kin
double epot = 0;			// Total potential energy
double virial = 0;			// Total virial


double dt = 0;				// Length of one timestep
double t = 0;				// Current time
int nt, nequil, nproduct;	// Number of timesteps
double tmax = 0;			// Total simulation time


TAverager avg_temp;			// Averager for the temperature
TAverager avg_epot;			// Averager for the total energy
TAverager avg_ekin;			// Averager for the total energy
TAverager avg_etot;			// Averager for the total energy
TAverager avg_vir;			// Averager for the total virial

int bincount, nbins;
double binwidth;
int* hist = NULL;			// Histogramm for radial distribution function

FILE* outHistogram;			// Outputfile for histogram
FILE* outTrajectories;		// Outputfile for trajectories
FILE* outAverages;			// Outputfile for averages


// int N, nequil, nproduct, naccepted, naccepted_prev, ntries, ntries_prev, bincount, nbins;
// double binwidth;


int main(int argc, char *argv[])
{

	if (argc < 9)
	{
		printf("Usage: %s <N> <eta> <dr> <nequil> <nproduct> <binwidth>\n", argv[0]);
		printf("\t<N>        Number of particles\n");
		printf("\t<rho>      Particle density\n");
		printf("\t<T>        Temperature\n");
		printf("\t<rc>       Lennard-Jones cut-off radius\n");
		printf("\t<dt>       Length of one timestep\n");
		printf("\t<nequil>   Number of equilibration timesteps to be performed\n");
		printf("\t<nproduct> Number of production timesteps to be performed\n");
		printf("\t<binwidth> Width of a bin for correlation length histogram\n");
		exit(EXIT_SUCCESS);
	}

	outHistogram = fopen("outHistogram.txt", "w+");
	outTrajectories = fopen("outTrajectories.txt", "w+");
	outAverages = fopen("outAverages.txt", "w+");

	// Parse commandline parameters
	N = atoi(argv[1]);
	rho = atof(argv[2]);
	T = atof(argv[3]);
	rc = atof(argv[4]);
	rc2 = rc*rc;
	dt = atof(argv[5]);
	nequil = atof(argv[6]);
	nproduct = atof(argv[7]);
	nt = nequil+nproduct;
	double tequil = nequil*dt;
	double tproduct = nproduct*dt;
	tmax = nt*dt;

	// Calculate corresponding system parameters
	lx = ly = sqrt((double)N/rho);

	printf("========== PARAMETERS ==========\n");
	printf("Particles:\t\t%d\n", N);
	printf("Density:\t\t%g\n", rho);
	printf("Simulationbox lx:\t%g\n", lx);
	printf("Simulationbox ly:\t%g\n", ly);
	printf("Temperature:\t\t%g\n", T);
	printf("Timestep length:\t%g\n", dt);
	printf("Equilibration steps:\t%d\n", nequil);
	printf("Production steps:\t%d\n", nproduct);
	printf("Total steps:\t%d\n", nt);
	printf("Equilibration time:\t%g\n", tequil);
	printf("Production time:\t%g\n", tproduct);
	printf("Total time:\t\t%g\n", tmax);
	printf("================================\n\n");


	printf("======== INIT PARTICLES ========\n");
	// Initialize arrays for particle positions & velocities
	r = new TVector[N];
	r_next = new TVector[N];
	v = new TVector[N];
	v_next = new TVector[N];

	// Put all particles equally spaced in the box and assign random velocities
	int nrows = sqrt(N);
	double dlx = lx/(double)nrows;
	double dly = ly/(double)nrows;
	vsum = .0;
	vsum2 = 0;
	for (int i = 0; i < N; i++)
	{
		// Positions
		r[i].x = i%nrows*dlx+0.5;
		r[i].y = floor(i/nrows)*dly+0.5;
		// Velocities
		v[i].x = rand_value(-1., 1.);
		v[i].y = rand_value(-1., 1.);
		vsum += v[i];
		vsum2 += v[i].x*v[i].x + v[i].y*v[i].y;
	}
	printf("Center of mass velocity after initialization:\t(%g,%g)\n", vsum.x, vsum.y);
	printf("Kinetic energy after initialization:\t\t%g\n", vsum2);
	printf("Instantaneous temperature after initialization:\t%g\n", vsum2/(2.0*(double)N));
	// Calculate average velocities
	vsum = vsum/(double)N;
	// Scalefactor for velocities to match the desired temperature (we neglect the fact
	// that the whole system with constrained center of mass has only (2N - 2) degrees
	// of freedom and approximate (2N - 2) \approx 2N since we won't run the simulation
	// with less than N = 100 particles.)
	double fs = sqrt(2.0*(double)N*T/vsum2);
	printf("Scaling factor for velocities:\t\t\t%g\n", fs);
	TVector vsumcheck;
	vsumcheck = .0;
	vsum2 = 0;
	for (int i = 0; i < N; i++)
	{
		v[i] = (v[i]-vsum)*fs;
		vsumcheck += v[i];
		vsum2 += v[i].x*v[i].x + v[i].y*v[i].y;
	}
	printf("Center of mass velocity after scaling:\t\t(%g,%g)\n", vsumcheck.x, vsumcheck.y);
	printf("Kinetic energy after scaling:\t\t\t%g\n", vsum2);
	printf("Instantaneous temperature after scaling:\t%g\n", vsum2/(2.0*(double)N));
	print_coords("outCoords_start.txt");
	printf("================================\n\n");


	printf("======== INIT POTENTIAL ========\n");
	// Init the potential
	init_lj_shift();
	F = new TVector[N];
	printf("Potential initialized.\n");
	printf("U(r_c)\t\t= %g\n", u_lj_shift);
	printf("U'(r_c)\t\t= %g\n", u_lj_deriv_shift);
	printf("U_s(r_c)\t= %g\n", u_lj_shifted(rc2));
	printf("U'_s(r_c)\t= %g\n", u_lj_deriv_shifted(rc2));
	printf("================================\n\n");


	printf("======== INIT AVERAGERS ========\n");
	avg_temp.init();
	avg_epot.init();
	avg_ekin.init();
	avg_etot.init();
	avg_vir.init();
	printf("Averagers initialized!\n");

	// Histogram for pair correlation function
	binwidth = atof(argv[8]);				// Width of a histogram-bin
	nbins = ceil(grmax/binwidth);			// Maximum correlation length to be measured should be L/2 due to periodic BC
	bincount = 0;							// Number of counts done on the histogram
	hist = new int[nbins];
	for (int i = 0; i < nbins; i++) {
		hist[i] = 0;
	}
	printf("Using histogram with %d bins of width %f\n", nbins, binwidth);
	printf("================================\n\n");

	printf("======= START INTEGRATION ======\n");
	t = 0;
	fprintf(outTrajectories, "t\tn\tr_x\t\tr_y\t\tv_x\t\tv_y\n");
	fprintf(outAverages, "t\tT(t)\t\t<T(t)>\t\tE_tot(T)\t\t<E_tot(T)>\n");
	for (int n = 0; n < nt; n++)
	{
		if (n == 0)
		{
			printf("Equilibration phase started.\n");
		}
		if (n == nequil)
		{
			printf("Production phase started.\n");
		}

		// Current time
		t = dt*n;
		if(debug) printf("t:\t%6.3f\t\n", t);
		// Calculate all forces
		forces();
		vsum = .0;
		vsum2 = .0;
		// update all particles
		for (int i = 0; i < N; i++)
		{
			// perform leap-frog-integration
			v_next[i] = v[i] + F[i]*dt;
			r_next[i] = r[i] + v_next[i]*dt;
			// Calculate energies
			vsum += v_next[i];
			// vsum2 += v[i].x*v[i].x + v[i].y*v[i].y; // naiv?
			vsum2 += pow(v_next[i].x+v[i].x, 2)/4.0 + pow(v_next[i].y+v[i].y, 2)/4.0; // sophisticated by Frenkel/Smit
			// update particle coordinates
			v[i] = v_next[i];
			r[i] = r_next[i];

			// Write trajectories to a file
			fprintf(outTrajectories, "%6.3f\t%6d\t%e\t%e\t%e\t%e\n", t, i, r[i].x, r[i].y, v[i].x, v[i].y);
		}

		// Equilibration phase, scale velocities to keep temperature
		if (n < nequil)
		{
			// Rescale velocities every ?? timesteps
			if (n%10 == 0)
			{
				scale_velocities();
			}
		}
		else if (n%nsamp == 0)
		{
			double Tt = vsum2/(2.0*(double)N);
			avg_temp.add(Tt);

			avg_epot.add(epot);
			avg_vir.add(virial);

			double ekin = 0.5*vsum2;
			avg_ekin.add(ekin);

			double etot = (epot + ekin);
			avg_etot.add(etot);

			update_histogram();

			fprintf(outAverages, "%6.3f\t%e\t%e\t%e\t%e\n", t, Tt, avg_temp.average(), etot, avg_etot.average());
		}

		if ((n+1)%(nt/10) == 0 || n == 0) {
			printf("Finished %5d (t = %5.1f) out of %d (t = %g) timesteps: %3.f %% <T> = %g\n", n+1, t, nt, tmax, (double)n/(double)nt*100, avg_temp.average());
		}

		if(debug) printf("\n");
	}
	printf("================================\n\n");

	print_coords("outCoords_end.txt");

	printf("Printing histogram for g(r)\n");
	for(int i = 0; i < nbins; i++) {
		double R = i*binwidth;
		double area = 2.0*PI*R*binwidth;
		// Multiply g(r) by two, since in the histogram we only counted each pair once, but each pair
		// gives two contributions to g(r)
		fprintf(outHistogram, "%f\t%f\n", R, 2.0*(double)hist[i]/(rho*area*(double)bincount*N));
	}

	delete [] r;
	delete [] r_next;
	delete [] v;
	delete [] v_next;
	delete [] F;
	r = r_next = v = v_next = F = NULL;

	// Close filepointer
	fclose(outHistogram); outHistogram = NULL;
	fclose(outTrajectories); outTrajectories = NULL;
	fclose(outAverages); outAverages = NULL;

	exit(EXIT_SUCCESS);
}

/** \brief Returns a uniformly distributed random number.
	
		Returns a uniformly distributed random number from the
		interval [min,max].
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param double Start of interval
	\param double End of interval
	\return double Random number
	\sa
**/
double rand_value(double min, double max) {
	return myran.doub() * (max - min) + min;
}

/** \brief Calculates distance between two particles
	
		Takes into account periodic boundary conditions and does not
		put particles back into the simulation box.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TVector First particle
	\param TVector Second particle
	\return TPartDist distance between particles
	\sa
**/
TPartDist distance_special(TVector vector1, TVector vector2)
{
	TPartDist pd;

	pd.dx = vector1.x - vector2.x;
	pd.dx = pd.dx - lx*round(pd.dx/lx);

	pd.dy = vector1.y - vector2.y;
	pd.dy = pd.dy - ly*round(pd.dy/ly);

	pd.r2 = pd.dx*pd.dx + pd.dy*pd.dy;

	return pd;
}

/** \brief Print out current coordinates and velocities to a file
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-12
	\param char Filename
	\return void
	\sa
**/
void print_coords(const char* filename) {
	FILE *outCoords;
	outCoords = fopen(filename, "w+");
	for (int i = 0; i < N; i++) {
		fprintf(outCoords, "%d\t%e\t%e\t%e\t%e\n", i, r[i].x, r[i].y, v[i].x, v[i].y);
	}
	fclose(outCoords); outCoords = NULL;
}


/** \brief Initialize shifted LJ-potential
	
		Calculates and stores the value of the LJ-potential at the cutoff
		radius.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\return void
	\sa
**/
void init_lj_shift()
{
	u_lj_shift = u_lj(rc2);
	u_lj_deriv_shift = u_lj_deriv(rc2);
}

/** \brief Calculate potential energy between two LJ-particles
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TPartDist Distance between two particles
	\return double Potential energy
	\sa
**/
double u_lj(double r2)
{
	double r2i = 1.0/r2;		// 1/r^2
	double r6i = r2i*r2i*r2i;	// 1/r^6
	
	return 4.0*r6i*(r6i-1.0);
}
/** \brief Calculate potential energy between two LJ-particles for shifted LJ-potential
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TPartDist Distance between two particles
	\return double Potential energy
	\sa
**/
double u_lj_shifted(double r2)
{
	return u_lj(r2) - u_lj_shift - (sqrt(r2) - rc)*u_lj_deriv_shift;
}


/** \brief Calculate the derivative of the LJ-potential
	
		If multiplied with a length, one gets the force.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TPartDist Distance between two particles
	\return double Derivative of the LJ-potential
	\sa
**/
double u_lj_deriv(double r2)
{
	double r2i = 1.0/r2;		// 1/r^2
	double r6i = r2i*r2i*r2i;	// 1/r^6

	return -48.0*r2i*r6i*(r6i-0.5);
}
/** \brief Calculate the derivative of the shifted LJ-potential
	
		If multiplied with a length, one gets the force.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TPartDist Distance between two particles
	\return double Derivative of the LJ-potential
	\sa
**/
double u_lj_deriv_shifted(double r2)
{
	return u_lj_deriv(r2) - u_lj_deriv_shift;
}

/** \brief Calculate the force between two LJ-particles with shifted LJ potential
	
		Calculates the force between two LJ-particles. The returned
		value containes the x- and y-component of the force.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\param TPartDist Distance between two particles
	\param double Virial between the two particles being considered
	\return TVector Force
	\sa
**/
TVector f_lj_shifted(TPartDist dr, double virij)
{
	TVector f;

	f.x = virij*dr.dx;
	f.y = virij*dr.dy;

	return f;
}

/** \brief Calculate the virial for two LJ-Particles
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-14
	\param TPartDist Distance between two particles
	\return double Virial for the two particles
	\sa
**/
double vir_lj_shifted(TPartDist dr)
{
	return -1.0*u_lj_deriv_shifted(dr.r2);
}

/** \brief Calculate the forces for all particle-pairs
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-13
	\return void
	\sa
**/
void forces()
{
	epot = 0.0;
	virial = 0.0;
	// Reset all forces
	for (int i = 0; i < N; i++)
	{
		F[i] = 0;
	}
	// Calculate current forces, loop over all pairs
	for (int i = 0; i < N-1; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			TPartDist rij = distance_special(r[i], r[j]);

			// Only evaluate if particles are within interaction distance
			if (rij.r2 > rc2)
			{
				continue;
			}

			double virij = vir_lj_shifted(rij);

			TVector flj = f_lj_shifted(rij, virij);
			F[i] += flj;
			F[j] -= flj;
			epot += u_lj_shifted(rij.r2);
			virial += virij;
		}
	}
}

/** \brief Rescale velocities
	
		Rescale velocities during equilibration phase to keep the temperature
		approximately fixed to the desired value.
	
	\author Christoph Tavan TU Berlin
	\date 2011-02-14
	\return void
	\sa
**/
void scale_velocities()
{
	double v2 = 0;
	for (int i = 0; i < N; i++)
	{
		v2 += v[i].x*v[i].x + v[i].y*v[i].y;
	}
	double fs = sqrt(2.0*(double)N*T/v2);
	for (int i = 0; i < N; i++)
	{
		v[i] *= fs;
	}
}


void update_histogram() {
	// Calculate histogram (correlation length)
	for (int i = 0; i < N-1; i++) {
		for (int j = i+1; j < N; j++) {
			// Get distance between particles i and j
			TPartDist rij = distance_special(r[i], r[j]);
			double r = sqrt(rij.r2);
			// Find correct bin of histogram
			int bin = floor(r/binwidth);
			if (bin >= nbins) {
				continue;
			}
			// printf("Bin:\t%d\n", bin);
			hist[bin]++;
		}
	}
	// Number of times averaging is done
	bincount++;
}