from igem_parts.bacterial_parts import j23101, j23101_seq, podd_backbone, podd_backbone_seq
import sbol3
import tyto
from typing import Tuple
from sbol_utilities.component import dna_component_with_sequence
from sbol_utilities.helper_functions import is_plasmid


def part_in_backbone(identity: str, part: sbol3.Component, backbone: sbol3.Component, linear:bool=False, **kwargs) -> Tuple[sbol3.Component, sbol3.Sequence]:
    """Creates a Part in Backbone Component and its Sequence.

    :param identity: The identity of the Component. The identity of Sequence is also identity with the suffix '_seq'.
    :param part: Part to be located in the backbone as SBOL Component.
    :param backbone: Backbone in wich the part is located as SBOL Component.
    :param linear: Boolean than indicates if the backbone is linear, by default it is seted to Flase which means that it has a circular topology.
    :param kwargs: Keyword arguments of any other Component attribute.
    :return: A tuple of Component and Sequence.
    """
    # check that backbone has a plasmid vector or child ontology term
    if is_plasmid(backbone)==False:
        raise TypeError('The backbone has no valid plasmid vector or child role')
    # check that the backbone and part has one sequence
    if len(backbone.sequences)!=1:
        raise ValueError(f'The backbone should have only one sequence, found {len(backbone.sequences)} sequences')
    if len(part.sequences)!=1:
        raise ValueError(f'The part should have only one sequence, found{len(part.sequences)} sequences')
    # check that the the last feature of backbone has 2 locations
    if len(backbone.features[-1].locations)!=2:
        raise ValueError(f'The backbone last feature should be the open backbone and should contain 2 Locations, found {len(backbone.features[-1].locations)} Locations')
    # get backbone sequence
    backbone_sequence = backbone.sequences[0].lookup().elements
    # compute open backbone sequences
    open_backbone_sequence_from_location1=backbone_sequence[backbone.features[-1].locations[0].start -1 : backbone.features[-1].locations[0].end]
    open_backbone_sequence_from_location2=backbone_sequence[backbone.features[-1].locations[1].start -1 : backbone.features[-1].locations[1].end]
    # extract part sequence
    part_sequence = part.sequences[0].lookup().elements
    # make new component sequence
    if linear:
        part_in_backbone_seq_str = open_backbone_sequence_from_location1 + part_sequence + open_backbone_sequence_from_location2
        topology_type = sbol3.SO_LINEAR
    else:
        part_in_backbone_seq_str = part_sequence + open_backbone_sequence_from_location2 + open_backbone_sequence_from_location1
        topology_type = sbol3.SO_CIRCULAR
    # part in backbone Component
    part_in_backbone_component, part_in_backbone_seq = dna_component_with_sequence(identity, part_in_backbone_seq_str, **kwargs)
    part_in_backbone_component.roles.append(tyto.SO.plasmid_vector) #review
    # defining Location
    part_subcomponent_location = sbol3.Range(sequence=part_in_backbone_seq, start=1, end=len(part_sequence))
    backbone_subcomponent_location = sbol3.Range(sequence=part_in_backbone_seq, start=len(part_sequence)+1, end=len(part_in_backbone_seq_str))
    source_location = sbol3.Range(sequence=backbone_sequence, start=backbone.features[-1].locations[0].start, end=backbone.features[-1].locations[0].end) # review
    # creating and adding features
    part_subcomponent = sbol3.SubComponent(part, roles=[tyto.SO.engineered_insert], locations=[part_subcomponent_location], role_integration='http://sbols.org/v3#mergeRoles')
    backbone_subcomponent = sbol3.SubComponent(backbone, locations=[backbone_subcomponent_location], source_locations=[source_location])  #[backbone.features[2].locations[0]]) #generalize source location
    part_in_backbone_component.features.append(part_subcomponent)
    part_in_backbone_component.features.append(backbone_subcomponent)
    # adding topology
    part_in_backbone_component.types.append(topology_type)
    return part_in_backbone_component, part_in_backbone_seq

doc = sbol3.Document()
sbol3.set_namespace('https://github.com/Gonza10V')
doc.add([j23101, j23101_seq, podd_backbone, podd_backbone_seq])
  
dist_j23101_in_podd, dist_j23101_in_podd_seq = part_in_backbone('j23101_in_podd', j23101, podd_backbone)

