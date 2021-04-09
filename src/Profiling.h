
struct Profile_Node;

struct Timed_Block_Signature
{
	const char* block_name;
	int line_number;
	const char* file_name;
};

struct Timed_Block
{	
	Profile_Node* profile_node;
	Timer timer;
	Timed_Block(const char* name, int line_number, const char* file_name);
	~Timed_Block();
};

struct Profile_Node
{
	Profile_Node* parent;
	int number_of_children;
	Profile_Node* children[16];
	Timed_Block_Signature signature;
	long unsigned int total_block_time;
};

struct Profile
{
	int number_of_nodes;
	int number_of_nodes_in_use;
	Profile_Node* nodes;
	Profile_Node* currently_executing_block;
};

#define TIMED_FUNCTION Timed_Block t_##__LINE__(__FUNCTION__, __LINE__, __FILE__)
#define TIMED_BLOCK(name) Timed_Block t_##__LINE__(name, __LINE__, __FILE__)

void init_profiling();
void print_profile();
