#include "catch.hpp"

#include "Attribute.hpp"
#include "CompareOperators.hpp"

TEST_CASE("cmp op test", "[CmpOps]") {

    Attribute::IsEqCmpOp isEq;

    REQUIRE(isEq(Attribute::Value(float(1)),Attribute::Value(float(1))));
    REQUIRE(!isEq(Attribute::Value(float(1.0002)),Attribute::Value(float(1.0001))));

    REQUIRE(isEq(Attribute::Value(int(1)),Attribute::Value(int(1))));
    REQUIRE(!isEq(Attribute::Value(int(1)),Attribute::Value(int(2))));

    Attribute::IsNotEqCmpOp isNotEq;

    REQUIRE(!isNotEq(Attribute::Value(float(1)),Attribute::Value(float(1))));
    REQUIRE(isNotEq(Attribute::Value(float(1.0002)),Attribute::Value(float(1.0001))));

    REQUIRE(!isNotEq(Attribute::Value(int(1)),Attribute::Value(int(1))));
    REQUIRE(isNotEq(Attribute::Value(int(1)),Attribute::Value(int(2))));

    Attribute::IsLessCmpOp isLess;

    REQUIRE(isLess(Attribute::Value(int(1)),Attribute::Value(int(2))));
    REQUIRE(!isLess(Attribute::Value(int(1)),Attribute::Value(int(1))));

    Attribute::IsLessOrEqCmpOp isLessOrEq;
    REQUIRE(isLessOrEq(Attribute::Value(int(1)),Attribute::Value(int(1))));

    Attribute::IsGreaterCmpOp isGreater;
    REQUIRE(isGreater(Attribute::Value(int(2)),Attribute::Value(int(1))));
    REQUIRE(!isGreater(Attribute::Value(int(2)),Attribute::Value(int(2))));

    Attribute::IsGreaterOrEqCmpOp isGreaterOrEq;
    REQUIRE(isGreaterOrEq(Attribute::Value(int(2)),Attribute::Value(int(2))));

    Attribute::Value val1(int(1));
    Attribute::Value val2(int(2));

    REQUIRE(val1 == val1);
    REQUIRE(val1 != val2);
}
