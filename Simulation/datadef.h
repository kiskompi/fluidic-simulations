enum Cell{
    C_F    =  0x0010,   /* This cell is a fluid cell */
    C_B    =  0x0000,   /* This cell is an obstacle/boundary cell */
    B_N    =  0x0001,   /* This obstacle cell has a fluid cell to the north */
    B_S    =  0x0002,   /* This obstacle cell has a fluid cell to the south */
    B_W    =  0x0004,   /* This obstacle cell has a fluid cell to the west */
    B_E    =  0x0008,   /* This obstacle cell has a fluid cell to the east */
    B_NW   =  (B_N | B_W),
    B_SW   =  (B_S | B_W),
    B_NE   =  (B_N | B_E),
    B_SE   =  (B_S | B_E),
    B_NSEW =  (B_N | B_S | B_E | B_W),
};


template<typename LT, typename RT>
constexpr inline unsigned int operator |(LT lval, RT rval) {
    return static_cast<unsigned int>(static_cast<unsigned int>(lval) | static_cast<unsigned int>(rval));
}

template<typename LT, typename RT>
constexpr inline unsigned int operator &(LT lval, RT rval) {
    return static_cast<unsigned int>(static_cast<unsigned int>(lval) & static_cast<unsigned int>(rval));
}

/* Macros for poisson(), denoting whether there is an obstacle cell
 * adjacent to some direction
 */
#define eps_E ((flag[i+1][j] & C_F)?1:0)
#define eps_W ((flag[i-1][j] & C_F)?1:0)
#define eps_N ((flag[i][j+1] & C_F)?1:0)
#define eps_S ((flag[i][j-1] & C_F)?1:0)
