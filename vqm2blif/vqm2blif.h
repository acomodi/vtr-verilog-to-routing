/*********************************************************************************************
*	vqm2blif.h
* This file contains useful structures and functions used by vqm2blif.cpp to
* parse a VQM netlist in order to output a corresponding BLIF netlist.
*
*			VQM to BLIF Convertor V.1.0
*
* Author:	S. Whitty
*		June 13, 2011
*
*********************************************************************************************/

#ifndef VQM2BLIF_H
#define VQM2BLIF_H

//============================================================================================
//				DEFINES
//============================================================================================
//preprocessor instruction to compile in debug mode
//outputs several intermediate files
#define TARGET_VPR

#define MAX_LEN 350
#define MAX_PPL 10
//============================================================================================
//				INCLUDES
//============================================================================================
#include <vector>
#include <string>
#include <map>

#include "../libvqm/vqm_dll.h"	//VQM Parser

extern "C" {
	#include "../libvpr/include/read_xml_arch_file.h"	//Architecture Parser
}

using namespace std;

//============================================================================================
//				GLOBAL TYPEDEFS & ENUMS
//============================================================================================

typedef vector<t_boolean> boolvec;
typedef vector<t_pin_def*> pinvec;
typedef vector<struct s_blif_subckt> scktvec;
typedef map <string,t_node_port_association*> portmap;
typedef pair<string,t_node_port_association*> portpair;

enum v_OptionBaseToken
{
	OT_VQM,
	OT_ARCH,
	OT_OUT,
	OT_DEBUG,
	OT_VERBOSE,
	OT_UNKNOWN
};

typedef map <string, v_OptionBaseToken> tokmap;
typedef pair <string, v_OptionBaseToken> tokpair;

//============================================================================================
//				DATA STRUCTURES
//============================================================================================
 
typedef struct s_blif_subckt{
	string inst_name;
	
	t_model* model_type;
	
	portmap input_cnxns;
	
	portmap output_cnxns;
	
} t_blif_subckt;
/* 				Subcircuit Structure
 *	-> used to represent a complete .subckt declaration of a .blif file.
 *
 *  inst_name:
 *	The instantiation name specified in the vqm file, preserved for readability 
 *	in the comments of the .blif
 *
 *  model_type:
 *	From VPR's architecture parser, contains data of the model's ports, 
 *	parameters, and modes. 
 *
 *  cnxns:
 *	Describes the intraconnect of the subckt. Maps a port from the model structure
 *	(logic_types.h, VPR) to an association containing pin data (vqm_dll.h, VQM). 
 *	If a port is left open, it will be mapped to the "hbpad" t_node_port_association.
 *	NOTE: Must map the port's name (string) because a port may be a bus, with multiple
 *	wires having the same portname. All assignments in BLIF must be 1-bit wide, so these
 *	buses must be flattened and mapped independantly.
 */
 
typedef struct s_blif_model{
	//General Variables
	string name;
	
	//Port Variables
	pinvec input_ports;
	
	pinvec output_ports;
	
	pinvec clock_ports; //NOTE: VPR currently doesn't support .clock, but could in the future
	
	//Assignment Variables
	int num_assignments;
	t_assign** array_of_assignments;
	
	//Subcircuit Variables
	scktvec subckts;
		
	//Logic (LUT) Variables
	//FUTURE WORK: Elaborate lut_masks into LUTs in the model
	
} t_blif_model;
/* 				Model Structure
 *	-> used to represent a complete .model declaration of a .blif file.
 *
 *  General Variables
 *  	name:
 * 	 	technical name of the model (e.g. ".model <name>").
 *
 *  Port Variables
 *	num_ports:
 *		total number of ports. 
 *	input_ports, output_ports:
 *		vectors containing references to the preallocated pins of the module. 
 *
 *  Assignment Variables
 *	num_assignments:
 *		number of assignment declarations in the model (e.g. "names ACCout~reg[0] ACCout[0] \n1 1").
 *	array_of_assignments:
 *		preallocated array containing assignment data, no need to reorganize.
 *
 *  Subcircuit Variables
 *	subckts:
 *		vector containing instances of t_blif_subckt structures and all necessary data.
 *
 *  Logic (LUT) Variables
 *	NOTE: 	Current version doesn't support soft logic declarations, only blackboxes. 
 *			data can be derived from lutmasks of lcells to allow downstream optimization. 
 *			Will investigate.
 */

#endif
