int32_t model_init(char * filename_arg, param_t * param_arg)
{
    // it is an error if both filename_arg and param_arg are supplies

    // XXX
    if (filename_arg) {
        // filename is supplied, 
        // read the model param and particle array from the file
        strcpy(filename, filename_arg);
        mode_read_file();
    } else {
        // init param and verify
        param = (param_arg ? *param_arg : param_default);
        verify_param();

        // allocate memory for particle array

        // init particle array filling the chamber with 
        // a distribution of particles, each particle has
        // a random position within the chamber and a random 
        // velocity direction, the velocity magnitude is set
        // to 300 degrees kelvin
    }

    // allocate and init location array

    // allocate and init radius array

    // register exit handler ??

}

int32_t model_read_file(void)
{
    // read file header and verify

    // read param from file

    // allocate memory for particle array

    // read particle array from file
}

int32_t model_write_file(void)
{
}

int32_t model_start(void)
{
}

int32_t model_stop(void)
{
}
