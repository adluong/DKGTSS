// SPDX-License-Identifier: LGPL-3.0-only
// This file is LGPL3 Licensed
pragma solidity ^0.8.0;

/**
 * @title Elliptic curve operations on twist points for alt_bn128
 * @author Mustafa Al-Bassam (mus@musalbas.com)
 * @dev Homepage: https://github.com/musalbas/solidity-BN256G2
 */

library BN256G2 {
    uint256 internal constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    uint256 internal constant TWISTBX = 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5;
    uint256 internal constant TWISTBY = 0x9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2;
    uint internal constant PTXX = 0;
    uint internal constant PTXY = 1;
    uint internal constant PTYX = 2;
    uint internal constant PTYY = 3;
    uint internal constant PTZX = 4;
    uint internal constant PTZY = 5;

    /**
     * @notice Add two twist points
     * @param pt1xx Coefficient 1 of x on point 1
     * @param pt1xy Coefficient 2 of x on point 1
     * @param pt1yx Coefficient 1 of y on point 1
     * @param pt1yy Coefficient 2 of y on point 1
     * @param pt2xx Coefficient 1 of x on point 2
     * @param pt2xy Coefficient 2 of x on point 2
     * @param pt2yx Coefficient 1 of y on point 2
     * @param pt2yy Coefficient 2 of y on point 2
     * @return (pt3xx, pt3xy, pt3yx, pt3yy)
     */
    function ECTwistAdd(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) public view returns (
        uint256, uint256,
        uint256, uint256
    ) {
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            if (!(
                pt2xx == 0 && pt2xy == 0 &&
                pt2yx == 0 && pt2yy == 0
            )) {
                assert(_isOnCurve(
                    pt2xx, pt2xy,
                    pt2yx, pt2yy
                ));
            }
            return (
                pt2xx, pt2xy,
                pt2yx, pt2yy
            );
        } else if (
            pt2xx == 0 && pt2xy == 0 &&
            pt2yx == 0 && pt2yy == 0
        ) {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
            return (
                pt1xx, pt1xy,
                pt1yx, pt1yy
            );
        }

        assert(_isOnCurve(
            pt1xx, pt1xy,
            pt1yx, pt1yy
        ));
        assert(_isOnCurve(
            pt2xx, pt2xy,
            pt2yx, pt2yy
        ));

        uint256[6] memory pt3 = _ECTwistAddJacobian(
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            1,     0,
            pt2xx, pt2xy,
            pt2yx, pt2yy,
            1,     0
        );

        return _fromJacobian(
            pt3[PTXX], pt3[PTXY],
            pt3[PTYX], pt3[PTYY],
            pt3[PTZX], pt3[PTZY]
        );
    }

    /**
     * @notice Multiply a twist point by a scalar
     * @param s     Scalar to multiply by
     * @param pt1xx Coefficient 1 of x
     * @param pt1xy Coefficient 2 of x
     * @param pt1yx Coefficient 1 of y
     * @param pt1yy Coefficient 2 of y
     * @return (pt2xx, pt2xy, pt2yx, pt2yy)
     */
    function ECTwistMul(
        uint256 s,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy
    ) public view returns (
        uint256, uint256,
        uint256, uint256
    ) {
        uint256 pt1zx = 1;
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            pt1xx = 1;
            pt1yx = 1;
            pt1zx = 0;
        } else {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
        }

        uint256[6] memory pt2 = _ECTwistMulJacobian(
            s,
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            pt1zx, 0
        );

        return _fromJacobian(
            pt2[PTXX], pt2[PTXY],
            pt2[PTYX], pt2[PTYY],
            pt2[PTZX], pt2[PTZY]
        );
    }

    /**
     * @notice Get the field modulus
     * @return The field modulus
     */
    function GetFieldModulus() public pure returns (uint256) {
        return FIELD_MODULUS;
    }

    function submod(uint256 a, uint256 b, uint256 n) internal pure returns (uint256) {
        return addmod(a, n - b, n);
    }

    function _FQ2Mul(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256, uint256) {
        return (
            submod(mulmod(xx, yx, FIELD_MODULUS), mulmod(xy, yy, FIELD_MODULUS), FIELD_MODULUS),
            addmod(mulmod(xx, yy, FIELD_MODULUS), mulmod(xy, yx, FIELD_MODULUS), FIELD_MODULUS)
        );
    }

    function _FQ2Muc(
        uint256 xx, uint256 xy,
        uint256 c
    ) internal pure returns (uint256, uint256) {
        return (
            mulmod(xx, c, FIELD_MODULUS),
            mulmod(xy, c, FIELD_MODULUS)
        );
    }

    function _FQ2Add(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256, uint256) {
        return (
            addmod(xx, yx, FIELD_MODULUS),
            addmod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Sub(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256 rx, uint256 ry) {
        return (
            submod(xx, yx, FIELD_MODULUS),
            submod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Div(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal view returns (uint256, uint256) {
        (yx, yy) = _FQ2Inv(yx, yy);
        return _FQ2Mul(xx, xy, yx, yy);
    }

    function _FQ2Inv(uint256 x, uint256 y) internal view returns (uint256, uint256) {
        uint256 inv = _modInv(addmod(mulmod(y, y, FIELD_MODULUS), mulmod(x, x, FIELD_MODULUS), FIELD_MODULUS), FIELD_MODULUS);
        return (
            mulmod(x, inv, FIELD_MODULUS),
            FIELD_MODULUS - mulmod(y, inv, FIELD_MODULUS)
        );
    }

    function _isOnCurve(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (bool) {
        uint256 yyx;
        uint256 yyy;
        uint256 xxxx;
        uint256 xxxy;
        (yyx, yyy) = _FQ2Mul(yx, yy, yx, yy);
        (xxxx, xxxy) = _FQ2Mul(xx, xy, xx, xy);
        (xxxx, xxxy) = _FQ2Mul(xxxx, xxxy, xx, xy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, xxxx, xxxy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, TWISTBX, TWISTBY);
        return yyx == 0 && yyy == 0;
    }

    function _modInv(uint256 a, uint256 n) internal view returns (uint256 result) {
        bool success;
        assembly {
            let freemem := mload(0x40)
            mstore(freemem, 0x20)
            mstore(add(freemem,0x20), 0x20)
            mstore(add(freemem,0x40), 0x20)
            mstore(add(freemem,0x60), a)
            mstore(add(freemem,0x80), sub(n, 2))
            mstore(add(freemem,0xA0), n)
            success := staticcall(sub(gas(), 2000), 5, freemem, 0xC0, freemem, 0x20)
            result := mload(freemem)
        }
        require(success);
    }

    function _fromJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal view returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) {
        uint256 invzx;
        uint256 invzy;
        (invzx, invzy) = _FQ2Inv(pt1zx, pt1zy);
        (pt2xx, pt2xy) = _FQ2Mul(pt1xx, pt1xy, invzx, invzy);
        (pt2yx, pt2yy) = _FQ2Mul(pt1yx, pt1yy, invzx, invzy);
    }

    function _ECTwistAddJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy) internal pure returns (uint256[6] memory pt3) {
            if (pt1zx == 0 && pt1zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt2xx, pt2xy,
                    pt2yx, pt2yy,
                    pt2zx, pt2zy
                );
                return pt3;
            } else if (pt2zx == 0 && pt2zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy
                );
                return pt3;
            }

            (pt2yx,     pt2yy)     = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // U1 = y2 * z1
            (pt3[PTYX], pt3[PTYY]) = _FQ2Mul(pt1yx, pt1yy, pt2zx, pt2zy); // U2 = y1 * z2
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // V1 = x2 * z1
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1xx, pt1xy, pt2zx, pt2zy); // V2 = x1 * z2

            if (pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY]) {
                if (pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY]) {
                    (
                        pt3[PTXX], pt3[PTXY],
                        pt3[PTYX], pt3[PTYY],
                        pt3[PTZX], pt3[PTZY]
                    ) = _ECTwistDoubleJacobian(pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy);
                    return pt3;
                }
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    1, 0,
                    1, 0,
                    0, 0
                );
                return pt3;
            }

            (pt2zx,     pt2zy)     = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // W = z1 * z2
            (pt1xx,     pt1xy)     = _FQ2Sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY]); // U = U1 - U2
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY]); // V = V1 - V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1yx, pt1yy, pt1yx,     pt1yy);     // V_squared = V * V
            (pt2yx,     pt2yy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY]); // V_squared_times_V2 = V_squared * V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1zx, pt1zy, pt1yx,     pt1yy);     // V_cubed = V * V_squared
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // newz = V_cubed * W
            (pt2xx,     pt2xy)     = _FQ2Mul(pt1xx, pt1xy, pt1xx,     pt1xy);     // U * U
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt2zx,     pt2zy);     // U * U * W
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt1zx,     pt1zy);     // U * U * W - V_cubed
            (pt2zx,     pt2zy)     = _FQ2Muc(pt2yx, pt2yy, 2);                    // 2 * V_squared_times_V2
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt2zx,     pt2zy);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
            (pt3[PTXX], pt3[PTXY]) = _FQ2Mul(pt1yx, pt1yy, pt2xx,     pt2xy);     // newx = V * A
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2yx, pt2yy, pt2xx,     pt2xy);     // V_squared_times_V2 - A
            (pt1yx,     pt1yy)     = _FQ2Mul(pt1xx, pt1xy, pt1yx,     pt1yy);     // U * (V_squared_times_V2 - A)
            (pt1xx,     pt1xy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY]); // V_cubed * U2
            (pt3[PTYX], pt3[PTYY]) = _FQ2Sub(pt1yx, pt1yy, pt1xx,     pt1xy);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
    }

    function _ECTwistDoubleJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy
    ) {
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 3);            // 3 * x
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1xx, pt1xy); // W = 3 * x * x
        (pt1zx, pt1zy) = _FQ2Mul(pt1yx, pt1yy, pt1zx, pt1zy); // S = y * z
        (pt2yx, pt2yy) = _FQ2Mul(pt1xx, pt1xy, pt1yx, pt1yy); // x * y
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // B = x * y * S
        (pt1xx, pt1xy) = _FQ2Mul(pt2xx, pt2xy, pt2xx, pt2xy); // W * W
        (pt2zx, pt2zy) = _FQ2Muc(pt2yx, pt2yy, 8);            // 8 * B
        (pt1xx, pt1xy) = _FQ2Sub(pt1xx, pt1xy, pt2zx, pt2zy); // H = W * W - 8 * B
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt1zx, pt1zy); // S_squared = S * S
        (pt2yx, pt2yy) = _FQ2Muc(pt2yx, pt2yy, 4);            // 4 * B
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt1xx, pt1xy); // 4 * B - H
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt2xx, pt2xy); // W * (4 * B - H)
        (pt2xx, pt2xy) = _FQ2Muc(pt1yx, pt1yy, 8);            // 8 * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1yx, pt1yy); // 8 * y * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt2zx, pt2zy); // 8 * y * y * S_squared
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt2xx, pt2xy); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 2);            // 2 * H
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // newx = 2 * H * S
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt2zx, pt2zy); // S * S_squared
        (pt2zx, pt2zy) = _FQ2Muc(pt2zx, pt2zy, 8);            // newz = 8 * S * S_squared
    }

    function _ECTwistMulJacobian(
        uint256 d,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (uint256[6] memory pt2) {
        while (d != 0) {
            if ((d & 1) != 0) {
                pt2 = _ECTwistAddJacobian(
                    pt2[PTXX], pt2[PTXY],
                    pt2[PTYX], pt2[PTYY],
                    pt2[PTZX], pt2[PTZY],
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy);
            }
            (
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            ) = _ECTwistDoubleJacobian(
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            );

            d = d / 2;
        }
    }
}
// This file is MIT Licensed.
//
// Copyright 2017 Christian Reitwiessner
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
pragma solidity ^0.8.0;
library Pairing {
    struct G1Point {
        uint X;
        uint Y;
    }
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }
    /// @return the generator of G1
    function P1() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }
    /// @return the generator of G2
    function P2() pure internal returns (G2Point memory) {
        return G2Point(
            [10857046999023057135944570762232829481370756359578518086990519993285655852781,
             11559732032986387107991004021392285783925812861821192530917403151452391805634],
            [8495653923123431417604973247489272438418190587263600148770280649306958101930,
             4082367875863433681332203403145435568316851327593401208105741076214120093531]
        );
    }
    /// @return the negation of p, i.e. p.addition(p.negate()) should be zero.
    function negate(G1Point memory p) pure internal returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p.X == 0 && p.Y == 0)
            return G1Point(0, 0);
        return G1Point(p.X, q - (p.Y % q));
    }
    /// @return r the sum of two points of G1
    function addition(G1Point memory p1, G1Point memory p2) internal view returns (G1Point memory r) {
        uint[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;
        assembly {
            success := staticcall(sub(gas(), 2000), 6, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
    }


    /// @return r the sum of two points of G2
    function addition(G2Point memory p1, G2Point memory p2) internal view returns (G2Point memory r) {
        (r.X[0], r.X[1], r.Y[0], r.Y[1]) = BN256G2.ECTwistAdd(p1.X[0],p1.X[1],p1.Y[0],p1.Y[1],p2.X[0],p2.X[1],p2.Y[0],p2.Y[1]);
    }


    /// @return r the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    function scalar_mul(G1Point memory p, uint s) internal view returns (G1Point memory r) {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := staticcall(sub(gas(), 2000), 7, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success);
    }
    /// @return the result of computing the pairing check
    /// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
    /// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
    /// return true.
    function pairing(G1Point[] memory p1, G2Point[] memory p2) internal view returns (bool) {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[1];
            input[i * 6 + 3] = p2[i].X[0];
            input[i * 6 + 4] = p2[i].Y[1];
            input[i * 6 + 5] = p2[i].Y[0];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := staticcall(sub(gas(), 2000), 8, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
        return out[0] != 0;
    }
    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2) internal view returns (bool) {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for three pairs.
    function pairingProd3(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2
    ) internal view returns (bool) {
        G1Point[] memory p1 = new G1Point[](3);
        G2Point[] memory p2 = new G2Point[](3);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for four pairs.
    function pairingProd4(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2,
            G1Point memory d1, G2Point memory d2
    ) internal view returns (bool) {
        G1Point[] memory p1 = new G1Point[](4);
        G2Point[] memory p2 = new G2Point[](4);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p1[3] = d1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        p2[3] = d2;
        return pairing(p1, p2);
    }
}

contract Verifier {
   using Pairing for *;
    struct VerifyingKey {
        Pairing.G2Point h;
        Pairing.G1Point g_alpha;
        Pairing.G2Point h_beta;
        Pairing.G1Point g_gamma;
        Pairing.G2Point h_gamma;
        Pairing.G1Point[] query;
    }
    struct Proof {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G1Point c;
    }
    //===========================================================================
    //Functions for paper
    struct Record{
        uint256[2]      z;
        uint256[2]      pk;
        uint256[2]      c;
        uint            stake;
        uint256         x;
        bool[2]         stt;
        bool[2]         rp;
    }

    mapping (address => Record)     l1; //list user
    address[]                       l6; //list address

    uint256[2][]                    pk;
    uint256[]                       sk;
    uint256[]                       x;
    uint256[3][]                    z;
    uint256[3][]                    hash;

    event ListZ(uint256[3][] list_z);
    event ListPK(uint256[2][] list_pkf);
    event ListShare(uint256[] list_share);
    event ListSK(uint256[] list_sk_frag);
    event ListHash(uint256[3][]);
    // event Rec(Record node);

    uint constant public stake_amount = 1 ether;
    uint public current_time = block.timestamp;
    uint public time1 = 1702467018;
    uint public time2 = 1702467018;
    
    function verifyC1(Proof memory proof, uint[4] memory input) payable public returns(bool){
        //require the caller has not called this function yet
        require(msg.value >= stake_amount);
        require(l1[msg.sender].stt[0] == false);
        uint[] memory inputValues = new uint[](4);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof, 1) == 0) {
            l1[msg.sender].z = [inputValues[0],inputValues[1]];
            l1[msg.sender].c = [inputValues[2],inputValues[3]];
            l1[msg.sender].stake = msg.value;
            //deposit is withdrawable by x
            l1[msg.sender].rp[1] = true;
            //done PInit
            l1[msg.sender].stt[0] = true;
            z.push([uint256(uint160(msg.sender)),inputValues[0],inputValues[1]]);
            return true;
        } else {
            return false;
        }
    }
    
    function verifyC2(Proof memory proof, uint[4] memory input) public returns(bool){
        require(l1[msg.sender].stt[0] == true);
        require(l1[msg.sender].stt[1] == false);
        uint[] memory inputValues = new uint[](4);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof, 2) == 0) {
            l1[msg.sender].pk = [inputValues[0],inputValues[1]];
            pk.push([inputValues[0],inputValues[1]]);
            hash.push([uint256(uint160(msg.sender)), inputValues[2],inputValues[3]]);
            l1[msg.sender].rp[0] = true;
            //DInit done
            l1[msg.sender].stt[1] = true;
            return true;
        } else {
            return false;
        }
    }//end skInit

    


    function report(address target, Proof memory proof, uint[3] memory input) public returns(bool){
        require((l1[target].pk[0] == input[1] && l1[target].pk[1] == input[2]) || (l1[target].c[0] == input[1] && l1[target].c[1] == input[2]));
        require(l1[target].rp[0] == true || l1[target].rp[1] == true);
        uint[] memory inputValues = new uint[](3);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof, 3) == 0) {
            l1[target].stake = 0;
            l1[msg.sender].stake += l1[target].stake/2;
            l1[msg.sender].rp = [false,false];
            l1[target].rp = [false,false];
            return true;
        } else {
            return false;
        }
    }

    //lai proof nua a?
    function submitSK(Proof memory proof, uint[3] memory input) public returns (bool){
        require(block.timestamp >= time1);
        require(l1[msg.sender].pk[0] == input[1]);  //pkf x
        require(l1[msg.sender].pk[1] == input[2]);  //pkf y
        uint[] memory inputValues = new uint[](3);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof,3) == 0) {
            sk.push(input[0]);
            l1[msg.sender].rp[0] == false;
            return true;
        } else {
            return false;
        }
    }

    function submitX(Proof memory proof, uint[3] memory input) public returns (bool){
        require(block.timestamp >= time2);
        require(l1[msg.sender].c[0] == input[1]);  //pkf x
        require(l1[msg.sender].c[1] == input[2]);  //pkf y
        uint[] memory inputValues = new uint[](3);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof,3) == 0) {
            // l1[msg.sender].x = inputValues[0];  //add share to the list
            x.push(input[0]);
            l1[msg.sender].rp[1] == false;
            return true;
        } else {
            return false;
        }
    }


    function query(uint option) public returns(bool){
        //Return list z. Called by dealers
        // uint256[][] memory buf;
        if(option == 1){
            emit ListZ(z);
            return true;
        }
        //Return list public keys. Called by all participants to reconstruct the MPK
        else if(option == 2){
            emit ListPK(pk);
            return true;
        }
        //Return secret keys and x-values. Called by at least k players.
         else if(option == 3){
            emit ListSK(sk);
            return true;
        }
        //Return secret keys and x-values. Called by at least k players.
         else if(option == 4){
            emit ListShare(x);
            return true;
        }
        else if(option == 5){
            emit ListHash(hash);
            return true;
        }
        return false;
    }

    //===========================================================================
    //
    //C1 pk
    function verifyingKey1() pure internal returns (VerifyingKey memory vk) {
        vk.h= Pairing.G2Point([uint256(0x00af8d3197d079ebac955c31555348a93648f435d1b182575e9b759f48dcafd0), uint256(0x0c0d4d856e1dc7da62702ad05749beafd2931c58399bab0d48239b8f0ca26f3b)], [uint256(0x20b65c18e7c35b380e282e8c000c0bc2a284113a9c3bedec7e9b72979d6fb6f7), uint256(0x16f84ec8d520e2f644498be7b490c5e8c42b1429fe0a3722c085f0a28b226068)]);
        vk.g_alpha = Pairing.G1Point(uint256(0x174127b2aafc85738e203791c58494faf1cc5efb5c5b7a09ddc9466ca1bdbf2f), uint256(0x27a3ff3d6a3e75fe48a87ee21b739775dded319eab43afc34ae3dd72b009c5eb));
        vk.h_beta = Pairing.G2Point([uint256(0x280ce8f096dbcc238707fd8cf9fe06da4a7adc2de3bebf9a0e747beaae2de318), uint256(0x1259d654e013f9c5b95d921f4c36d59cdc80174ff8900e7f689d35c19fa50db1)], [uint256(0x16bdb160479ea18a4511d872fbf5ba1e628fa8593a7678d427106fb939f77304), uint256(0x2e255b4ea3482a0bb33836007f41bb66404601a5ffa2c2ca745cf9d082cdb9cc)]);
        vk.g_gamma = Pairing.G1Point(uint256(0x1b6e5de41fb909d147e1227e20cf3f11ee1fede20107173c8a691ac996c6f10b), uint256(0x06078eb05e1555855b6fb07c904d4553acc731db277c34057eee1b05e1ea650e));
        vk.h_gamma = Pairing.G2Point([uint256(0x00af8d3197d079ebac955c31555348a93648f435d1b182575e9b759f48dcafd0), uint256(0x0c0d4d856e1dc7da62702ad05749beafd2931c58399bab0d48239b8f0ca26f3b)], [uint256(0x20b65c18e7c35b380e282e8c000c0bc2a284113a9c3bedec7e9b72979d6fb6f7), uint256(0x16f84ec8d520e2f644498be7b490c5e8c42b1429fe0a3722c085f0a28b226068)]);
        vk.query = new Pairing.G1Point[](5);
        vk.query[0] = Pairing.G1Point(uint256(0x090eacbc53d8bce66830042be7b63683c2f7d352916e9add5e545c4bccc4cb51), uint256(0x2a9362876e9d5fc0172026ab9ca369b50d42a89593bb4f47f18d79d35d18b240));
        vk.query[1] = Pairing.G1Point(uint256(0x29ec848bee6325ef7ba23e9f8cad289115d524e511725cc894c46a5cd570acf0), uint256(0x25f70fbbddb416aade9389a2b1c4b99f75a0d86b643b3f950c1443a9184d5567));
        vk.query[2] = Pairing.G1Point(uint256(0x19f7b1d19d0e9b875b8773ab4786beaf238cfaa0a2efc117e89e831f10bd6554), uint256(0x0704834961ac79eb6d42cca15c71d4a934447732810bd59f686f4cc0cc198286));
        vk.query[3] = Pairing.G1Point(uint256(0x11e3a81f10aae2d85ce6292672b5d4475f8086b9e44cc431baff297cf2308fb6), uint256(0x0b0c12e274eadc8a7705074f1a26d8ce4006344de6bf25efbe9209d21f7ad1cb));
        vk.query[4] = Pairing.G1Point(uint256(0x0273dd09422d3c88e84497ae4b9640aad26bcc44efbbd2322c1d20db860c2134), uint256(0x26c863c2698e66f434e6d1b50a05356367b3c6fe235a280232b7ec695c828903));
    }

    //C2 pk
    function verifyingKey2() pure internal returns (VerifyingKey memory vk) {
        vk.h= Pairing.G2Point([uint256(0x15bc968aaa9123547ea6bc46d471cc5208137287aa8e2b2e01db48bc435a480d), uint256(0x2f0b2d47bc3ee39c85cc88ede75af7e725c437ff3ddb3aee26e81ed4fb32d27c)], [uint256(0x02db2dbb66168a5b0f0566ea1c17648b4afd929152ff245d82ef5234110bccf5), uint256(0x2a759ccabbb426430467423aef6d69648a5c1c3afa4b62e6e238acc854874139)]);
        vk.g_alpha = Pairing.G1Point(uint256(0x2af47df2e08ed6bb83b479641fa054228bfc895d1ebb53a4a51eec4020c7c2e0), uint256(0x22ca4df22a0d31ca27794c06c2612b2ba0c4dca80fb75cc8d569cca10d17b7f0));
        vk.h_beta = Pairing.G2Point([uint256(0x05d86c2c8f23db56258d840b9b232d1adb14548ed0ca538ba75255367b45eee5), uint256(0x0073c398bc9d5140dfdcf327d51ba3523b2ff98c5117c52a4b0c242eebe509b3)], [uint256(0x2fb450d0d86d493dcacc3900aa9f90c2f8d7e5f0722b183cf65e944c5a3c8615), uint256(0x25552b592c18272b3f192b561e0058381f034e7b7d9895e16ba12d0e15937441)]);
        vk.g_gamma = Pairing.G1Point(uint256(0x1315538addcb8094d9a713c947b7512d5891299993cb778d66e72e954fd8b61c), uint256(0x1cf28dee920d0623f34ee31d5c2a64fd7b95f29929ba01a2dfe71529bfc6141f));
        vk.h_gamma = Pairing.G2Point([uint256(0x15bc968aaa9123547ea6bc46d471cc5208137287aa8e2b2e01db48bc435a480d), uint256(0x2f0b2d47bc3ee39c85cc88ede75af7e725c437ff3ddb3aee26e81ed4fb32d27c)], [uint256(0x02db2dbb66168a5b0f0566ea1c17648b4afd929152ff245d82ef5234110bccf5), uint256(0x2a759ccabbb426430467423aef6d69648a5c1c3afa4b62e6e238acc854874139)]);
        vk.query = new Pairing.G1Point[](5);
        vk.query[0] = Pairing.G1Point(uint256(0x1bea7736256691654005ec39dbe4a8030697877894d6b8589d3ca1688c43db1d), uint256(0x1a5330cfc6b70d9282b56af8356e037561f3b5b122def4aeff1631961b2a059b));
        vk.query[1] = Pairing.G1Point(uint256(0x15204ee535c9e94e3be0dd8150ed5da25c6aca0d899763e1445485d9b6ae8a8c), uint256(0x2dc1cc210cacac5508337b7f6f573b14c97b167df0179cc4256a3ae15c6f403c));
        vk.query[2] = Pairing.G1Point(uint256(0x27604e1155b1fc79efee1cdd5f0c109b0e3b5a9447e8ca127338b634ba48d4a8), uint256(0x121ae9206e2b9873cc5895aefa06d138e788850858d7c1bc7dcd8f5c2195668f));
        vk.query[3] = Pairing.G1Point(uint256(0x198191e5a94442c976448fa8758defd287ccf1c266ca74948705cb71d8208039), uint256(0x1d85ee643deda1e58d1569e008acdddc8c467966ad0564922183e9ee78ea4a1a));
        vk.query[4] = Pairing.G1Point(uint256(0x1c68cfc6ca1569553b9aaebb74f9a08b4b3b1ec93d38aafb5e7adde5b63c2f48), uint256(0x140c964d44cb27abd011dfd6fbb73384f92059571893c5d7bbcc0e664876bda2));
    }

    function verifyingKey3() pure internal returns (VerifyingKey memory vk) {
        vk.h= Pairing.G2Point([uint256(0x260dea6724029f2fa4a9a8e9a25801d2c73e6b5c35cfb16ed27bc699069ea2ab), uint256(0x11ecf7c05d3d2d58244481184bdf8e5fc7c97705e1e4c6326d3eb9d3d222a601)], [uint256(0x03abbf2588774fc741951fa698e862714f2eb1a2f9bafa3cb424586997b9510d), uint256(0x1677f892978ef8691f1ba45e24e262af48d1e762e82d7e4b1d18470c160b3ba6)]);
        vk.g_alpha = Pairing.G1Point(uint256(0x148b91c5737297dbf037c5a6b44427f3dfa3c1e3c746ed880bc207655ec715e5), uint256(0x2c3fa763e1f8ee200be41c269b79a76906e759c300490dbc1bc7d2d47487fba8));
        vk.h_beta = Pairing.G2Point([uint256(0x00a6097ddbdec13f8b94aa4e172dec7d4e1d6800e2430b62c6a22f159660ecde), uint256(0x1dd5e59d19fce78f40d216ec5d64f32cdda90da6ea7b1ef96a07e9bc3c3ff5e3)], [uint256(0x0430c32b7a94f054eeac0a8c68c5fd288e0be59043c82ca4af3759496c2d9b00), uint256(0x07c0d409819dc9cb6b44b49d9bc70a60adba28fdf722cd002f46dd27fad1602b)]);
        vk.g_gamma = Pairing.G1Point(uint256(0x1c41008d02ff84087fa292364010389bb4e3dd50b76c976cf6cff09b06eea12c), uint256(0x0b80d612ee3143621f0d4f9fc553a813756eb05f756e9e2706a3352cd160e2cc));
        vk.h_gamma = Pairing.G2Point([uint256(0x260dea6724029f2fa4a9a8e9a25801d2c73e6b5c35cfb16ed27bc699069ea2ab), uint256(0x11ecf7c05d3d2d58244481184bdf8e5fc7c97705e1e4c6326d3eb9d3d222a601)], [uint256(0x03abbf2588774fc741951fa698e862714f2eb1a2f9bafa3cb424586997b9510d), uint256(0x1677f892978ef8691f1ba45e24e262af48d1e762e82d7e4b1d18470c160b3ba6)]);
        vk.query = new Pairing.G1Point[](4);
        vk.query[0] = Pairing.G1Point(uint256(0x2e5911ceeab94b46e0a8dfbca62d34115d77defee3c8d973f0f8e39e52879f93), uint256(0x2b70ffd3fda6ba68dda1efdd00c3b78e285b61fa6b4a32c36349d1a983455129));
        vk.query[1] = Pairing.G1Point(uint256(0x04beaf7831e935318f39ccd866e7791ad5214618209788c085b125f7a74a426a), uint256(0x22fb370e12cac3f837beb2c4689d1954e5a666f326b6385e8021b4d030768302));
        vk.query[2] = Pairing.G1Point(uint256(0x18ed3127e832ebaaf362bad9f1acfa3810a698343efd3d934057b0b98455a7a2), uint256(0x07ace69e5191a820c1d16ac01209bfe840340fb4f10ed7586b06b56d00c5bf6d));
        vk.query[3] = Pairing.G1Point(uint256(0x048e245cf59d50a0f317b9f1623410909790fd4704246023f1928bbbaeffde83), uint256(0x152a30c79e3422dbf8f281fe75fafa038a4fe1b527d11d349f2a5b6f328019c2));
    }
    function verify(uint[] memory input, Proof memory proof, uint16 option) internal view returns (uint) {
        uint256 snark_scalar_field = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        VerifyingKey memory vk;
        if(option == 1){
            vk = verifyingKey1();
        }
        else if(option == 2){
            vk = verifyingKey2();
        }
        else{
            vk = verifyingKey3();
        }
        // VerifyingKey memory vk = (option == true) ? verifyingKey(): verifyingKey2();
        require(input.length + 1 == vk.query.length);
        // Compute the linear combination vk_x
        Pairing.G1Point memory vk_x = Pairing.G1Point(0, 0);
        for (uint i = 0; i < input.length; i++) {
            require(input[i] < snark_scalar_field);
            vk_x = Pairing.addition(vk_x, Pairing.scalar_mul(vk.query[i + 1], input[i]));
        }
        vk_x = Pairing.addition(vk_x, vk.query[0]);
        /**
         * e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma})
         *                              * e(C, H)
         * where psi = \sum_{i=0}^l input_i pvk.query[i]
         */
        if (!Pairing.pairingProd4(vk.g_alpha, vk.h_beta, vk_x, vk.h_gamma, proof.c, vk.h, Pairing.negate(Pairing.addition(proof.a, vk.g_alpha)), Pairing.addition(proof.b, vk.h_beta))) return 1;
        /**
         * e(A, H^{gamma}) = e(G^{gamma}, B)
         */
        if (!Pairing.pairingProd2(proof.a, vk.h_gamma, Pairing.negate(vk.g_gamma), proof.b)) return 2;
        return 0;
    }
}
