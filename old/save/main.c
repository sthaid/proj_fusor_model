common.h
- model file struct
- params
- global vars
    modelfile_t mf;

mf_max

model_run flag


// -----------------------------------------------------------

int32_t main(int32_t argc, char **argv)
{
    // initialize
    initialize();

    // display and control the model
    display_and_control();
}

int32_t initialize()
{
    // get options

    // init model
}

void display_and_control()
{
    // init sdl

    while (true) {
        // display params
        // display graphs
        // display chamber

        // read and process events
        // - list the events
    }
}

// --- seperate file ---

int32_t model_init()
{
    // if this is a new file then create the initial state
    // if modelfile doesnt exist then 
    // create and init empty modelfile
    // memory map the modelfile
    // validate modelfile 

    // set modelfile globals
    // - addr

    // register model_exit handler

    // create model work and control threads

}

void exit_handler()
{
    // pause the model
}

void * control_thread(void * cx)
{
}

void * work_thread(void * cx)
{
}

