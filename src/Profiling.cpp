Profile raytrace_profile;

bool operator==(Timed_Block_Signature s, Timed_Block_Signature t)
{
	return strcmp(s.block_name, t.block_name) == 0 && strcmp(s.file_name, t.file_name) == 0 && s.line_number == t.line_number;
}

void init_profiling()
{
	raytrace_profile = {};
	raytrace_profile.nodes = (Profile_Node*)alloc(256 * sizeof(Profile_Node));
	raytrace_profile.nodes[0].parent = NULL;
	for(int i = 0; i < 16; ++i) raytrace_profile.nodes[0].children[i] = NULL;
	raytrace_profile.nodes[0].signature = {"Raytrace profile", __LINE__, __FILE__};
	raytrace_profile.nodes[0].total_block_time = 0L;

	++raytrace_profile.number_of_nodes_in_use;
	
	raytrace_profile.currently_executing_block = raytrace_profile.nodes;
}

void print_profile_node(Profile_Node* node, int depth)
{

	long unsigned int node_time = node->total_block_time;
	for(int i = 0; i < node->number_of_children; ++i) node_time -= node->children[i]->total_block_time;
	printf("%*s| %s: %f s with %d calls\n", depth, " ", node->signature.block_name, cycles_to_s(node_time), node->number_of_calls);

	for(int i = 0; i < node->number_of_children; ++i) print_profile_node(node->children[i], depth + 1);
}

void print_profile()
{
	//printf("\t| %s: %lu cycles\n", node->signature.block_name, node->total_block_time);
	for(int i = 0; i < raytrace_profile.nodes[0].number_of_children; ++i)
	{
		print_profile_node(raytrace_profile.nodes[0].children[i], 1);
	}
}

bool profile_node_exists(Timed_Block_Signature signature)
{
	Profile_Node* currently_executing_block = raytrace_profile.currently_executing_block;
	for(int i = 0; i < currently_executing_block->number_of_children; ++i)
	{
		if(currently_executing_block->children[i]->signature == signature) return true;
	}
	return false;
}

Profile_Node* get_profile_node(Timed_Block_Signature signature)
{
	Profile_Node* currently_executing_block = raytrace_profile.currently_executing_block;
	for(int i = 0; i < currently_executing_block->number_of_children; ++i)
	{
		if(currently_executing_block->children[i]->signature == signature) return currently_executing_block->children[i];
	}
	return NULL;
}

Profile_Node* get_free_profile_node()
{
	Profile_Node* free_node = raytrace_profile.nodes + raytrace_profile.number_of_nodes_in_use;
	++raytrace_profile.number_of_nodes_in_use;
	return free_node;
}

void set_currently_executing_block(Profile_Node* block)
{
	raytrace_profile.currently_executing_block = block;
}

void insert_profile_node(Profile_Node* node)
{
	Profile_Node** new_child_node;
	for(	new_child_node = raytrace_profile.currently_executing_block->children;
		*new_child_node != NULL;
		++new_child_node);
	*new_child_node = node;
	++raytrace_profile.currently_executing_block->number_of_children;
}

Timed_Block::Timed_Block(const char* block_name, int line_number, const char* file_name)
{
	Timed_Block_Signature signature = {block_name, line_number, file_name};
	if(profile_node_exists(signature))
	{
		this->profile_node = get_profile_node(signature);
	}
	else
	{
		this->profile_node = get_free_profile_node();
		this->profile_node->signature = signature;
		insert_profile_node(this->profile_node);
	}
	this->profile_node->parent = raytrace_profile.currently_executing_block;

	set_currently_executing_block(this->profile_node);
	++this->profile_node->number_of_calls;
	start_timer(&this->timer);
}

Timed_Block::~Timed_Block()
{
	stop_timer(&this->timer);
	this->profile_node->total_block_time += elapsed_time_in_cycles(&this->timer);
	set_currently_executing_block(this->profile_node->parent);
}
