package element

import "github.com/shopspring/decimal"

type Element struct {
	AtomicNumber    int     // Atomic number of the element
	Symbol          string  // Symbol of the element (e.g., "H", "O", "He")
	Name            string  // Full name of the element (e.g., "Hydrogen", "Oxygen")
	AtomicWeight    decimal.Decimal // Atomic weight (usually a floating point number)
	Electronegativity float64 // Electronegativity (Pauling scale)
	VanDerWaalsRadius float64 // Van der Waals radius (in picometers, pm)
	Group           int     // The group number(Column) in the periodic table (1-18)
	Period          int     // The period number(Row) in the periodic table (1-7)
}

type PeriodicTable struct {
	Elements []Element
}

func NewPeriodicTable() *PeriodicTable {
    elements := []Element{
        {1, "H", "Hydrogen", decimal.NewFromFloat(1.008), 2.20, 120.0, 1, 1},
		{2, "He", "Helium", decimal.NewFromFloat(4.002602), 0.0, 140.0, 18, 1},
		{3, "Li", "Lithium", decimal.NewFromFloat(6.94), 0.98, 182.0, 1, 2},
		{4, "Be", "Beryllium", decimal.NewFromFloat(9.0122), 1.57, 153.0, 2, 2},
		{5, "B", "Boron", decimal.NewFromFloat(10.81), 2.04, 192.0, 13, 2},
		{6, "C", "Carbon", decimal.NewFromFloat(12.011), 2.55, 170.0, 14, 2},
		{7, "N", "Nitrogen", decimal.NewFromFloat(14.007), 3.04, 155.0, 15, 2},
		{8, "O", "Oxygen", decimal.NewFromFloat(15.999), 3.44, 152.0, 16, 2},
		{9, "F", "Fluorine", decimal.NewFromFloat(18.998403163), 3.98, 147.0, 17, 2},
		{10, "Ne", "Neon", decimal.NewFromFloat(20.1797), 0.0, 154.0, 18, 2},
		{11, "Na", "Sodium", decimal.NewFromFloat(22.98976928), 0.93, 180.0, 1, 3},
		{12, "Mg", "Magnesium", decimal.NewFromFloat(24.305), 1.31, 173.0, 2, 3},
		{13, "Al", "Aluminum", decimal.NewFromFloat(26.9815385), 1.61, 184.0, 13, 3},
		{14, "Si", "Silicon", decimal.NewFromFloat(28.085), 1.90, 210.0, 14, 3},
		{15, "P", "Phosphorus", decimal.NewFromFloat(30.973761998), 2.19, 175.0, 15, 3},
		{16, "S", "Sulfur", decimal.NewFromFloat(32.065), 2.58, 180.0, 16, 3},
		{17, "Cl", "Chlorine", decimal.NewFromFloat(35.45), 3.16, 175.0, 17, 3},
		{18, "Ar", "Argon", decimal.NewFromFloat(39.948), 0.0, 188.0, 18, 3},
		{19, "K", "Potassium", decimal.NewFromFloat(39.0983), 0.82, 275.0, 1, 4},
		{20, "Ca", "Calcium", decimal.NewFromFloat(40.078), 1.00, 231.0, 2, 4},
		{21, "Sc", "Scandium", decimal.NewFromFloat(44.955908), 1.36, 184.0, 3, 4},
		{22, "Ti", "Titanium", decimal.NewFromFloat(47.867), 1.54, 176.0, 4, 4},
		{23, "V", "Vanadium", decimal.NewFromFloat(50.9415), 1.63, 171.0, 5, 4},
		{24, "Cr", "Chromium", decimal.NewFromFloat(52.0), 1.66, 139.0, 6, 4},
		{25, "Mn", "Manganese", decimal.NewFromFloat(54.938044), 1.55, 161.0, 7, 4},
		{26, "Fe", "Iron", decimal.NewFromFloat(55.845), 1.83, 155.0, 8, 4},
		{27, "Co", "Cobalt", decimal.NewFromFloat(58.933194), 1.88, 152.0, 9, 4},
		{28, "Ni", "Nickel", decimal.NewFromFloat(58.6934), 1.91, 149.0, 10, 4},
		{29, "Cu", "Copper", decimal.NewFromFloat(63.546), 1.90, 135.0, 11, 4},
		{30, "Zn", "Zinc", decimal.NewFromFloat(65.38), 1.65, 139.0, 12, 4},
		{31, "Ga", "Gallium", decimal.NewFromFloat(69.723), 1.81, 187.0, 13, 4},
		{32, "Ge", "Germanium", decimal.NewFromFloat(72.63), 2.01, 211.0, 14, 4},
		{33, "As", "Arsenic", decimal.NewFromFloat(74.921595), 2.18, 185.0, 15, 4},
		{34, "Se", "Selenium", decimal.NewFromFloat(78.971), 2.55, 190.0, 16, 4},
		{35, "Br", "Bromine", decimal.NewFromFloat(79.904), 2.96, 185.0, 17, 4},
		{36, "Kr", "Krypton", decimal.NewFromFloat(83.798), 3.00, 202.0, 18, 4},
		{37, "Rb", "Rubidium", decimal.NewFromFloat(85.4678), 0.82, 303.0, 1, 5},
		{38, "Sr", "Strontium", decimal.NewFromFloat(87.62), 0.95, 249.0, 2, 5},
		{39, "Y", "Yttrium", decimal.NewFromFloat(88.90584), 1.22, 253.0, 3, 5},
		{40, "Zr", "Zirconium", decimal.NewFromFloat(91.224), 1.33, 200.0, 4, 5},
		{41, "Nb", "Niobium", decimal.NewFromFloat(92.90637), 1.60, 198.0, 5, 5},
		{42, "Mo", "Molybdenum", decimal.NewFromFloat(95.95), 2.16, 200.0, 6, 5},
		{43, "Tc", "Technetium", decimal.NewFromFloat(98), 2.00, 217.0, 7, 5},
		{44, "Ru", "Ruthenium", decimal.NewFromFloat(101.07), 2.20, 207.0, 8, 5},
		{45, "Rh", "Rhodium", decimal.NewFromFloat(102.90550), 2.28, 198.0, 9, 5},
		{46, "Pd", "Palladium", decimal.NewFromFloat(106.42), 2.20, 163.0, 10, 5},
		{47, "Ag", "Silver", decimal.NewFromFloat(107.8682), 1.93, 172.0, 11, 5},
		{48, "Cd", "Cadmium", decimal.NewFromFloat(112.411), 1.69, 158.0, 12, 5},
		{49, "In", "Indium", decimal.NewFromFloat(114.818), 1.78, 193.0, 13, 5},
		{50, "Sn", "Tin", decimal.NewFromFloat(118.710), 1.96, 217.0, 14, 5},
		{51, "Sb", "Antimony", decimal.NewFromFloat(121.760), 2.05, 202.0, 15, 5},
		{52, "Te", "Tellurium", decimal.NewFromFloat(127.60), 2.01, 206.0, 16, 5},
		{53, "I", "Iodine", decimal.NewFromFloat(126.90447), 2.66, 198.0, 17, 5},
		{54, "Xe", "Xenon", decimal.NewFromFloat(131.293), 2.60, 216.0, 18, 5},
		{55, "Cs", "Cesium", decimal.NewFromFloat(132.90545196), 0.79, 343.0, 1, 6},
		{56, "Ba", "Barium", decimal.NewFromFloat(137.327), 0.89, 253.0, 2, 6},
		{57, "La", "Lanthanum", decimal.NewFromFloat(138.90547), 1.10, 262.0, 3, 6},
		{58, "Ce", "Cerium", decimal.NewFromFloat(140.116), 1.12, 266.0, 3, 6},
		{59, "Pr", "Praseodymium", decimal.NewFromFloat(140.90766), 1.13, 267.0, 3, 6},
		{60, "Nd", "Neodymium", decimal.NewFromFloat(144.242), 1.14, 270.0, 3, 6},
		{61, "Pm", "Promethium", decimal.NewFromFloat(145), 1.13, 271.0, 3, 6},
		{62, "Sm", "Samarium", decimal.NewFromFloat(150.36), 1.17, 274.0, 3, 6},
		{63, "Eu", "Europium", decimal.NewFromFloat(151.964), 1.20, 277.0, 3, 6},
		{64, "Gd", "Gadolinium", decimal.NewFromFloat(157.25), 1.20, 280.0, 3, 6},
		{65, "Tb", "Terbium", decimal.NewFromFloat(158.92535), 1.23, 282.0, 3, 6},
		{66, "Dy", "Dysprosium", decimal.NewFromFloat(162.500), 1.22, 285.0, 3, 6},
		{67, "Ho", "Holmium", decimal.NewFromFloat(164.93033), 1.23, 287.0, 3, 6},
		{68, "Er", "Erbium", decimal.NewFromFloat(167.259), 1.24, 289.0, 3, 6},
		{69, "Tm", "Thulium", decimal.NewFromFloat(168.93422), 1.25, 292.0, 3, 6},
		{70, "Yb", "Ytterbium", decimal.NewFromFloat(173.04), 1.10, 294.0, 3, 6},
		{71, "Lu", "Lutetium", decimal.NewFromFloat(174.9668), 1.27, 296.0, 3, 6},
		{72, "Hf", "Hafnium", decimal.NewFromFloat(178.49), 1.30, 208.0, 4, 6},
		{73, "Ta", "Tantalum", decimal.NewFromFloat(180.94788), 1.50, 200.0, 5, 6},
		{74, "W", "Tungsten", decimal.NewFromFloat(183.84), 2.36, 193.0, 6, 6},
		{75, "Re", "Rhenium", decimal.NewFromFloat(186.207), 1.90, 188.0, 7, 6},
		{76, "Os", "Osmium", decimal.NewFromFloat(190.23), 2.20, 190.0, 8, 6},
		{77, "Ir", "Iridium", decimal.NewFromFloat(192.217), 2.20, 180.0, 9, 6},
		{78, "Pt", "Platinum", decimal.NewFromFloat(195.084), 2.28, 177.0, 10, 6},
		{79, "Au", "Gold", decimal.NewFromFloat(196.966569), 2.54, 144.0, 11, 6},
		{80, "Hg", "Mercury", decimal.NewFromFloat(200.592), 2.00, 155.0, 12, 6},
		{81, "Tl", "Thallium", decimal.NewFromFloat(204.38), 1.62, 196.0, 13, 6},
		{82, "Pb", "Lead", decimal.NewFromFloat(207.2), 2.33, 202.0, 14, 6},
		{83, "Bi", "Bismuth", decimal.NewFromFloat(208.98040), 2.02, 207.0, 15, 6},
		{84, "Po", "Polonium", decimal.NewFromFloat(209), 2.00, 202.0, 16, 6},
		{85, "At", "Astatine", decimal.NewFromFloat(210), 2.2, 202.0, 17, 6},
		{86, "Rn", "Radon", decimal.NewFromFloat(222), 2.2, 220.0, 18, 6},
		{87, "Fr", "Francium", decimal.NewFromFloat(223), 0.7, 330.0, 1, 7},
		{88, "Ra", "Radium", decimal.NewFromFloat(226), 0.9, 215.0, 2, 7},
		{89, "Ac", "Actinium", decimal.NewFromFloat(227), 1.1, 216.0, 3, 7},
		{90, "Th", "Thorium", decimal.NewFromFloat(232.03805), 1.3, 232.0, 3, 7},
		{91, "Pa", "Protactinium", decimal.NewFromFloat(231.03588), 1.5, 231.0, 4, 7},
		{92, "U", "Uranium", decimal.NewFromFloat(238.02891), 1.38, 244.0, 5, 7},
		{93, "Np", "Neptunium", decimal.NewFromFloat(237), 1.36, 259.0, 6, 7},
		{94, "Pu", "Plutonium", decimal.NewFromFloat(244), 1.28, 263.0, 7, 7},
		{95, "Am", "Americium", decimal.NewFromFloat(243), 1.13, 267.0, 8, 7},
		{96, "Cm", "Curium", decimal.NewFromFloat(247), 1.3, 273.0, 9, 7},
		{97, "Bk", "Berkelium", decimal.NewFromFloat(247), 1.3, 276.0, 10, 7},
		{98, "Cf", "Californium", decimal.NewFromFloat(251), 1.3, 281.0, 11, 7},
		{99, "Es", "Einsteinium", decimal.NewFromFloat(252), 1.5, 282.0, 12, 7},
		{100, "Fm", "Fermium", decimal.NewFromFloat(257), 1.6, 287.0, 13, 7},
		{101, "Md", "Mendelevium", decimal.NewFromFloat(258), 1.7, 290.0, 14, 7},
		{102, "No", "Nobelium", decimal.NewFromFloat(259), 1.7, 292.0, 15, 7},
		{103, "Lr", "Lawrencium", decimal.NewFromFloat(262), 1.7, 294.0, 16, 7},
		{104, "Rf", "Rutherfordium", decimal.NewFromFloat(267), 1.6, 297.0, 4, 7},
		{105, "Db", "Dubnium", decimal.NewFromFloat(270), 1.6, 300.0, 5, 7},
		{106, "Sg", "Seaborgium", decimal.NewFromFloat(271), 1.6, 303.0, 6, 7},
		{107, "Bh", "Bohrium", decimal.NewFromFloat(270), 1.6, 305.0, 7, 7},
		{108, "Hs", "Hassium", decimal.NewFromFloat(277), 1.6, 310.0, 8, 7},
		{109, "Mt", "Meitnerium", decimal.NewFromFloat(276), 1.6, 315.0, 9, 7},
		{110, "Ds", "Darmstadtium", decimal.NewFromFloat(281), 1.6, 318.0, 10, 7},
		{111, "Rg", "Roentgenium", decimal.NewFromFloat(280), 1.6, 320.0, 11, 7},
		{112, "Cn", "Copernicium", decimal.NewFromFloat(285), 1.6, 325.0, 12, 7},
		{113, "Nh", "Nihonium", decimal.NewFromFloat(284), 1.6, 330.0, 13, 7},
		{114, "Fl", "Flerovium", decimal.NewFromFloat(289), 1.6, 335.0, 14, 7},
		{115, "Mc", "Moscovium", decimal.NewFromFloat(288), 1.6, 340.0, 15, 7},
		{116, "Lv", "Livermorium", decimal.NewFromFloat(293), 1.6, 345.0, 16, 7},
		{117, "Ts", "Tennessine", decimal.NewFromFloat(294), 1.6, 350.0, 17, 7},
		{118, "Og", "Oganesson", decimal.NewFromFloat(294), 2.0, 360.0, 18, 7},
    }
    return &PeriodicTable{Elements:elements}
}

func (pt *PeriodicTable) FindElementBySymbol(symbol string) (*Element, bool) {
	for _, elem := range pt.Elements {
		if elem.Symbol == symbol {
			return &elem, true
		}
	}
	return nil, false //if the passed symbol isn't in the correct format it won't be found 
}


