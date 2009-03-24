#include <cstdio>
#include <cstdlib>
#include <symbols.hxx>

//const size_t offset(5321735);
const size_t offset(4915199);

int main(int argc, char ** argv) {

	if(argc < 2) {
		fprintf(stderr,
			"Usage: %s <restart file>\n", argv[0]);
		exit(1);
	} // if

	// open restart file
	FILE * restart_file = fopen(argv[1], "r+");
	if(restart_file == NULL) {
		fprintf(stderr, "Error opening %s\n", argv[1]);
		exit(1);
	} // if

	// seek to beginning of function pointers
	fseek(restart_file, offset, SEEK_SET);

	// write new function pointers
	fwrite(&new_field, sizeof(void *), 1, restart_file);
	fwrite(&delete_field, sizeof(void *), 1, restart_file);
	fwrite(&new_material_coefficients, sizeof(void *), 1, restart_file);
	fwrite(&delete_material_coefficients, sizeof(void *), 1, restart_file);
	fwrite(&advance_b, sizeof(void *), 1, restart_file);
	fwrite(&advance_e, sizeof(void *), 1, restart_file);
	fwrite(&energy_f, sizeof(void *), 1, restart_file);
	fwrite(&clear_jf, sizeof(void *), 1, restart_file);
	fwrite(&synchronize_jf, sizeof(void *), 1, restart_file);
	fwrite(&clear_rhof, sizeof(void *), 1, restart_file);
	fwrite(&synchronize_rho, sizeof(void *), 1, restart_file);
	fwrite(&compute_rhob, sizeof(void *), 1, restart_file);
	fwrite(&compute_curl_b, sizeof(void *), 1, restart_file);
	fwrite(&synchronize_tang_e_norm_b, sizeof(void *), 1, restart_file);
	fwrite(&compute_div_e_err, sizeof(void *), 1, restart_file);
	fwrite(&compute_rms_div_e_err, sizeof(void *), 1, restart_file);
	fwrite(&clean_div_e, sizeof(void *), 1, restart_file);
	fwrite(&compute_div_b_err, sizeof(void *), 1, restart_file);
	fwrite(&compute_rms_div_b_err, sizeof(void *), 1, restart_file);
	fwrite(&clean_div_b, sizeof(void *), 1, restart_file);

	// close files
	fclose(restart_file);

	return 0;
} // main
