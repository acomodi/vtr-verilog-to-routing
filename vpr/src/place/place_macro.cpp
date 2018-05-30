#include <cstdio>
#include <ctime>
#include <cmath>
#include <sstream>
#include <map>
using namespace std;

#include "vtr_assert.h"
#include "vtr_memory.h"
#include "vtr_util.h"

#include "vpr_types.h"
#include "vpr_error.h"
#include "physical_types.h"
#include "globals.h"
#include "place.h"
#include "read_xml_arch_file.h"
#include "place_macro.h"
#include "vpr_utils.h"
#include "echo_files.h"


/******************** File-scope variables declarations **********************/

/* f_idirect_from_blk_pin array allow us to quickly find pins that could be in a    *
 * direct connection. Values stored is the index of the possible direct connection  *
 * as specified in the arch file, OPEN (-1) is stored for pins that could not be    *
 * part of a direct chain conneciton.                                               *
 * [0...device_ctx.num_block_types-1][0...num_pins-1]                               */
static std::vector<std::vector<int>> f_idirect_from_blk_pin;

/* f_direct_type_from_blk_pin array stores the value SOURCE if the pin is the       *
 * from_pin, SINK if the pin is the to_pin in the direct connection as specified in *
 * the arch file, OPEN (-1) is stored for pins that could not be part of a direct   *
 * chain conneciton.                                                                *
 * [0...device_ctx.num_block_types-1][0...num_pins-1]                               */
static std::vector<std::vector<int>> f_direct_type_from_blk_pin;

/* f_imacro_from_blk_pin maps a blk_num to the corresponding macro index.           *
 * If the block is not part of a macro, the value OPEN (-1) is stored.              *
 * [0...cluster_ctx.clb_nlist.blocks().size()-1]                                    */
static vtr::vector_map<ClusterBlockId, int> f_imacro_from_iblk;


/******************** Subroutine declarations ********************************/

static void find_all_the_macro (int * num_of_macro, std::vector<ClusterBlockId> &pl_macro_member_blk_num_of_this_blk,
		std::vector<int> &pl_macro_idirect, std::vector<int> &pl_macro_num_members, std::vector<std::vector<ClusterBlockId>> &pl_macro_member_blk_num);

static void alloc_and_load_imacro_from_iblk(t_pl_macro * macros, int num_macros);

static void write_place_macros(std::string filename, const t_pl_macro* macros, int num_macros);

static bool is_constant_clb_net(ClusterNetId clb_net);

static bool net_is_driven_by_direct(ClusterNetId clb_net);

static void validate_macros(t_pl_macro* macros, int num_macro);
/******************** Subroutine definitions *********************************/


static void find_all_the_macro (int * num_of_macro, std::vector<ClusterBlockId> &pl_macro_member_blk_num_of_this_blk,
		std::vector<int> &pl_macro_idirect, std::vector<int> &pl_macro_num_members, std::vector<std::vector<ClusterBlockId>> &pl_macro_member_blk_num) {

	/* Compute required size:                                                *
	 * Go through all the pins with possible direct connections in           *
	 * f_idirect_from_blk_pin. Count the number of heads (which is the same  *
	 * as the number macros) and also the length of each macro               *
	 * Head - blocks with to_pin OPEN and from_pin connected                 *
	 * Tail - blocks with to_pin connected and from_pin OPEN                 */

    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& clb_nlist = cluster_ctx.clb_nlist;

	int num_macro = 0;
	for (auto blk_id : cluster_ctx.clb_nlist.blocks()) {
        vtr::printf("Checking Block %s for macro head\n", cluster_ctx.clb_nlist.block_name(blk_id).c_str());

        for (ClusterPinId sink_pin : clb_nlist.block_input_pins(blk_id)) {

            auto sink_itype = clb_nlist.block_type(blk_id)->index;
            auto sink_ipin = clb_nlist.pin_physical_index(sink_pin);
			ClusterNetId sink_net = clb_nlist.block_net(blk_id, sink_ipin);

			int sink_idirect = f_idirect_from_blk_pin[sink_itype][sink_ipin];

            vtr::printf("  Checking input pin: %s\n", clb_nlist.pin_name(sink_pin).c_str(), (sink_net) ? clb_nlist.net_name(sink_net).c_str() : "N/A");
			// Identify potential macro head blocks (i.e. start of a macro)
            //
            // The input SINK (to_pin) of a potential HEAD macro must have:
            //  * an incomming direct connection, and
            //  * either no connected net, or a net not driven by a direct connection
            //
            // Note that the restriction that the net is not driven from another direct ensures that
            // blocks in the middle of a chain with internal connections are not detected has potential
            // head blocks.

            if (sink_idirect == OPEN) continue; //Not a direct connection
            if (sink_net != ClusterNetId::INVALID() && net_is_driven_by_direct(sink_net)) continue; //Connected net via to direct

            for (ClusterPinId driver_pin : clb_nlist.block_output_pins(blk_id)) {

                auto driver_itype = clb_nlist.block_type(blk_id)->index;
                auto driver_ipin = clb_nlist.pin_physical_index(driver_pin);
                ClusterNetId driver_net = clb_nlist.block_net(blk_id, driver_ipin);

                int driver_idirect = f_idirect_from_blk_pin[driver_itype][driver_ipin];

                // Confirm whether this is a head macro
                //
                // The output SOURCE (from_pin) of a true head macro will:
                //  * drive another block with the same direct connection
                if (sink_idirect != driver_idirect) continue; //Different direct
                if (driver_net == ClusterNetId::INVALID()) continue; //No connected downstream blocks
                    
                vtr::printf("    Checking output pin: %s net: %s\n", clb_nlist.pin_name(driver_pin).c_str(), clb_nlist.net_name(driver_net).c_str());

                // Mark down that this is the first block in the macro
                pl_macro_member_blk_num_of_this_blk[0] = blk_id;
                pl_macro_idirect[num_macro] = sink_idirect;
                
                // Increment the num_member count.
                pl_macro_num_members[num_macro]++;
                
                // Also find out how many members are in the macros, 
                // there are at least 2 members - 1 head and 1 tail.
                
                // Initialize the variables
                ClusterNetId next_net_id = driver_net;
                ClusterBlockId next_blk_id = blk_id;

                // Start finding the other members
                while (next_net_id != ClusterNetId::INVALID()) {
                    ClusterNetId curr_net_id = next_net_id;
                    
                    // Assume that carry chains only has 1 sink - direct connection
                    VTR_ASSERT(cluster_ctx.clb_nlist.net_sinks(curr_net_id).size() == 1);
                    next_blk_id = cluster_ctx.clb_nlist.net_pin_block(curr_net_id, 1);			

                    if (cluster_ctx.clb_nlist.block_type(next_blk_id) != cluster_ctx.clb_nlist.block_type(blk_id)) {
                        //Different block types, end of chain
                        //
                        //Currently we detect only macros connected to the same block type
                        break;
                    }
                    
                    // Assume that the from_iblk_pin index is the same for the next block
                    VTR_ASSERT(f_idirect_from_blk_pin[cluster_ctx.clb_nlist.block_type(next_blk_id)->index][driver_ipin] == driver_idirect
                            && f_direct_type_from_blk_pin[cluster_ctx.clb_nlist.block_type(next_blk_id)->index][driver_ipin] == SOURCE);
                    next_net_id = cluster_ctx.clb_nlist.block_net(next_blk_id, driver_ipin);

                    // Mark down this block as a member of the macro
                    int imember = pl_macro_num_members[num_macro];
                    pl_macro_member_blk_num_of_this_blk[imember] = next_blk_id;

                    // Increment the num_member count.
                    pl_macro_num_members[num_macro]++;

                } // Found all the members of this macro at this point

                // Allocate the second dimension of the blk_num array since I now know the size
                pl_macro_member_blk_num[num_macro].resize(pl_macro_num_members[num_macro]);
                // Copy the data from the temporary array to the newly allocated array.
                for (int imember = 0; imember < pl_macro_num_members[num_macro]; imember++) {
                    pl_macro_member_blk_num[num_macro][imember] = pl_macro_member_blk_num_of_this_blk[imember];
                }

                // Increment the macro count
                num_macro++;
            } // Finish going through all the pins for from_pins.
		} // Finish going through all the pins for to_pins.
	} // Finish going through all blocks.
	
	// Now, all the data is readily stored in the temporary data structures.
	*num_of_macro = num_macro;
}


int alloc_and_load_placement_macros(t_direct_inf* directs, int num_directs, t_pl_macro ** macros){
	
	/* This function allocates and loads the macros placement macros   *
	 * and returns the total number of macros in 2 steps.              *
	 *   1) Allocate temporary data structure for maximum possible     *
	 *      size and loops through all the blocks storing the data     *
	 *      relevant to the carry chains. At the same time, also count *
	 *      the amount of memory required for the actual variables.    *
	 *   2) Allocate the actual variables with the exact amount of     *
	 *      memory. Then loads the data from the temporary data        *
	 *       structures before freeing them.                           *
	 *                                                                 *
	 * For pl_macro_member_blk_num, allocate for the first dimension   *
	 * only at first. Allocate for the second dimemsion when I know    *
	 * the size. Otherwise, the array is going to be of size           *
	 * cluster_ctx.clb_nlist.blocks().size()^2 (There are big		   *
	 * benckmarks VPR that have cluster_ctx.clb_nlist.blocks().size()  *
	 * in the 100k's range).										   *
	 *																   *
	 * The placement macro array is freed by the caller(s).            */

	/* Declaration of local variables */
	int imacro, imember, num_macro;
	auto& cluster_ctx = g_vpr_ctx.clustering();

	/* Allocate maximum memory for temporary variables. */
	std::vector<int> pl_macro_idirect(cluster_ctx.clb_nlist.blocks().size());
	std::vector<int> pl_macro_num_members(cluster_ctx.clb_nlist.blocks().size());
	std::vector<std::vector<ClusterBlockId>> pl_macro_member_blk_num(cluster_ctx.clb_nlist.blocks().size());
	std::vector<ClusterBlockId> pl_macro_member_blk_num_of_this_blk(cluster_ctx.clb_nlist.blocks().size());
	
	t_pl_macro * macro = nullptr;
	
	/* Sets up the required variables. */
	alloc_and_load_idirect_from_blk_pin(directs, num_directs, 
			f_idirect_from_blk_pin, f_direct_type_from_blk_pin);


	/* Compute required size:                                                *
	 * Go through all the pins with possible direct connections in           *
	 * f_idirect_from_blk_pin. Count the number of heads (which is the same  *
	 * as the number macros) and also the length of each macro               *
	 * Head - blocks with to_pin OPEN and from_pin connected                 *
	 * Tail - blocks with to_pin connected and from_pin OPEN                 */
	num_macro = 0;
	find_all_the_macro (&num_macro, pl_macro_member_blk_num_of_this_blk, 
			pl_macro_idirect, pl_macro_num_members, pl_macro_member_blk_num);

	/* Allocate the memories for the macro. */
	macro = (t_pl_macro *) vtr::malloc (num_macro * sizeof(t_pl_macro));

	/* Allocate the memories for the chaim members.             *
	 * Load the values from the temporary data structures.      */
	for (imacro = 0; imacro < num_macro; imacro++) {
		macro[imacro].num_blocks = pl_macro_num_members[imacro];
		macro[imacro].members = (t_pl_macro_member *) vtr::malloc(macro[imacro].num_blocks * sizeof(t_pl_macro_member));

		/* Load the values for each member of the macro */
		for (imember = 0; imember < macro[imacro].num_blocks; imember++) {
			macro[imacro].members[imember].x_offset = imember * directs[pl_macro_idirect[imacro]].x_offset;
			macro[imacro].members[imember].y_offset = imember * directs[pl_macro_idirect[imacro]].y_offset;
			macro[imacro].members[imember].z_offset = directs[pl_macro_idirect[imacro]].z_offset;
			macro[imacro].members[imember].blk_index = pl_macro_member_blk_num[imacro][imember];
		}
	}
	
	/* Returns the pointer to the macro by reference. */
	*macros = macro;

    if(isEchoFileEnabled(E_ECHO_PLACE_MACROS)) {
        write_place_macros(getEchoFileName(E_ECHO_PLACE_MACROS), *macros, num_macro);
    }

    validate_macros(*macros, num_macro);

	return (num_macro);
}

void get_imacro_from_iblk(int *imacro, ClusterBlockId iblk, t_pl_macro *macros, int num_macros) {

	/* This mapping is needed for fast lookup's whether the block with index *
	 * iblk belongs to a placement macro or not.                             *
	 *                                                                       *
	 * The array f_imacro_from_iblk is used for the mapping for speed reason *
	 * [0...cluster_ctx.clb_nlist.blocks().size()-1]                                                    */

	/* If the array is not allocated and loaded, allocate it.                */ 
	if (f_imacro_from_iblk.size() == 0) {
		alloc_and_load_imacro_from_iblk(macros, num_macros);
	}

	/* Return the imacro for the block. */
	*imacro = f_imacro_from_iblk[iblk];

}

/* Allocates and loads imacro_from_iblk array. */
static void alloc_and_load_imacro_from_iblk(t_pl_macro *macros, int num_macros) {
	int imacro, imember;
    auto& cluster_ctx = g_vpr_ctx.clustering();

	f_imacro_from_iblk.resize(cluster_ctx.clb_nlist.blocks().size());

	/* Allocate and initialize the values to OPEN (-1). */
	for (auto blk_id : cluster_ctx.clb_nlist.blocks()) {
		f_imacro_from_iblk.insert(blk_id, OPEN);
	}
	
	/* Load the values */
	for (imacro = 0; imacro < num_macros; imacro++) {
		for (imember = 0; imember < macros[imacro].num_blocks; imember++) {
			ClusterBlockId blk_id = macros[imacro].members[imember].blk_index;
			f_imacro_from_iblk.insert(blk_id, imacro);
		}
	}
}

void free_placement_macros_structs() {

	/* This function frees up all the static data structures used. */

	// This frees up the two arrays
    f_idirect_from_blk_pin.clear();
    f_direct_type_from_blk_pin.clear();
}

static void write_place_macros(std::string filename, const t_pl_macro *macros, int num_macros) {

    FILE* f = vtr::fopen(filename.c_str(), "w");

    fprintf(f, "#Identified Placement macros\n");
    fprintf(f, "Num_Macros: %d\n", num_macros);
    for (int imacro = 0; imacro < num_macros; ++imacro) {
        const t_pl_macro* macro = &macros[imacro];
        fprintf(f, "Macro_Id: %d, Num_Blocks: %d\n", imacro, macro->num_blocks);
        fprintf(f, "------------------------------------------------------\n");
        for (int imember = 0; imember < macro->num_blocks; ++imember) {
            const t_pl_macro_member* macro_memb = &macro->members[imember];
            fprintf(f, "Block_Id: %zu, x_offset: %d, y_offset: %d, z_offset: %d\n",
                    size_t(macro_memb->blk_index),
                    macro_memb->x_offset,
                    macro_memb->y_offset,
                    macro_memb->z_offset);
        }
        fprintf(f, "\n");
    }

    fprintf(f, "\n");

    fprintf(f, "#Macro-related direct connections\n");
    fprintf(f, "type      type_pin  is_direct direct_type\n");
    fprintf(f, "------------------------------------------\n");
    auto& device_ctx = g_vpr_ctx.device();
    for (int itype = 0; itype < device_ctx.num_block_types; ++itype) {
        t_type_descriptor* type = &device_ctx.block_types[itype];

        for (int ipin = 0; ipin < type->num_pins; ++ipin) {
            if (f_idirect_from_blk_pin[itype][ipin] != OPEN) {
                if (f_direct_type_from_blk_pin[itype][ipin] == SOURCE) {
                    fprintf(f, "%-9s %-9d true      SOURCE    \n", type->name, ipin);
                } else {
                    VTR_ASSERT(f_direct_type_from_blk_pin[itype][ipin] == SINK);
                    fprintf(f, "%-9s %-9d true      SINK      \n", type->name, ipin);
                }
            } else {
                VTR_ASSERT(f_direct_type_from_blk_pin[itype][ipin] == OPEN);
            }
        }

    }

    fclose(f);
}

static bool is_constant_clb_net(ClusterNetId clb_net) {
    auto& atom_ctx = g_vpr_ctx.atom();
    AtomNetId atom_net = atom_ctx.lookup.atom_net(clb_net);

    return atom_ctx.nlist.net_is_constant(atom_net);
}

static bool net_is_driven_by_direct(ClusterNetId clb_net) {
    auto& cluster_ctx = g_vpr_ctx.clustering();
    
	ClusterBlockId block_id = cluster_ctx.clb_nlist.net_driver_block(clb_net);
	int pin_index = cluster_ctx.clb_nlist.net_pin_physical_index(clb_net, 0);

    auto direct = f_idirect_from_blk_pin[cluster_ctx.clb_nlist.block_type(block_id)->index][pin_index];

    return direct != OPEN;
}

static void validate_macros(t_pl_macro* macros, int num_macros) {
    //Perform sanity checks on macros
    auto& cluster_ctx = g_vpr_ctx.clustering();

    //Verify that blocks only appear in a single macro
    std::multimap<ClusterBlockId,int> block_to_macro;
    for (int imacro = 0; imacro < num_macros; ++imacro) {
        for (int imember = 0; imember < macros[imacro].num_blocks; ++imember) {
            ClusterBlockId iblk = macros[imacro].members[imember].blk_index;

            block_to_macro.emplace(iblk, imacro);
        }
    }

    for (auto blk_id : cluster_ctx.clb_nlist.blocks()) {
        auto range = block_to_macro.equal_range(blk_id);

        int blk_macro_cnt = std::distance(range.first, range.second);
        if (blk_macro_cnt > 1) {
            std::stringstream msg;
            msg << "Block #" << size_t(blk_id) << " '" << cluster_ctx.clb_nlist.block_name(blk_id) << "'"
                << " appears in " << blk_macro_cnt << " placement macros (should appear in at most one). Related Macros:\n";

            for (auto iter = range.first; iter != range.second; ++iter) {
                int imacro = iter->second;
                msg << "  Macro #: " << imacro << "\n";
            }

            VPR_THROW(VPR_ERROR_PLACE, msg.str().c_str());
        }
    }
}
