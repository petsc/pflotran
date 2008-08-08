extern "C"
{
void assign_c_ptr_(void** return_arg, void* pointer_arg)
{
   *return_arg = pointer_arg;
}

void assign_c_array_ptr_(void** return_arg, void* pointer_arg)
{
   *return_arg = pointer_arg;
}

}
