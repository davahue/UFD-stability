# UFD-stability
This repository includes codes allowing the stability determination for different categories of Urban Flood Drifters (UFDs). It also includes the processed results for easier determination of the inception of movement.

The UFD subcategories considered are those of Table 2 of: Bayon, Valero and Franca (2023):
+------------------------------------+-------------------------+----+-----------------------+--------------------------------------------------------------------------------+
|                                    | Category                | ID | Subcategory           | Description                                                                    |
+------------------------------------+-------------------------+----+-----------------------+--------------------------------------------------------------------------------+
|                                    |                         | V1 | Two-wheelers          | Bikes, motorbikes and e-scooters.                                              |
|                                    | Vehicles                | V2 | Cars                  | Cars and other light four-wheel vehicles designed to transport of passengers.  |
| Typified UFDs                      | (UFD-V)                 | V3 | Vans                  | Vans and other heavy four-wheel vehicles designed to transport materials       |
|                                    |                         |    |                       | and stock.                                                                     |
|                                    |                         | V4 | Caravans & RVs        | Vehicles designed to provide habitable space (RV: recreational vehicle).       |
|                                    |                         | V5 | Large heavy vehicles  | Vehicles designed to transport a large amount of people or goods               |
|                                    |                         |    |                       | (buses, trucks, trains, boats, etc.).                                          |
|                                    |                         | F1 | Urban fixtures        | Facilities designed to provide a public service in streets                     |
|                                    | Furniture               |    |                       | (bins, waste containers, etc.).                                                |
|                                    | (UFD-F)                 | F2 | Household equipment   | Facilities from private front (and back) gardens that can be carried           |
|                                    |                         |    |                       | by floods (tanks, garden sheds, etc.).                                         |
+------------------------------------+-------------------------+----+-----------------------+--------------------------------------------------------------------------------+
|                                    |                         | DC | Construction          | Debris that can be dragged from construction sites or damaged buildings.       |
|                                    |                         | DM | Metal                 | Metal debris, predominantly of constructive origin (sheets, pipes, etc.).      |
| Heterog. UFDs                      | (UFD-H)                 | DP | Plastic               | Plastics and textile objects of small dimensions and irregular shape.          |
|                                    |                         | DW | Wood                  | Natural wood (trunks, branches, etc.) and processed wood.                      |
|                                    |                         | DO | Others                | Other drifters of uncertain origin (food, tableware, leaves, sediment, etc.).  |
+------------------------------------+-------------------------+----+-----------------------+--------------------------------------------------------------------------------+

- Bayón, Arnau, Valero, Daniel, & Franca,  Mário J. (2023). Urban Flood Drifters (UFDs): identification, classification and characterisation. ArXiv preprint: https://arxiv.org/abs/2304.01780
