
/*
CompareOperators.hpp - This file is part of BayesTyper (v0.9)


The MIT License (MIT)

Copyright (c) 2016 Jonas Andreas Sibbesen and Lasse Maretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#ifndef COMPARE_OPERATOR
#define COMPARE_OPERATOR

#include "boost/variant.hpp"

#include "Attribute.hpp"
#include "Utils.hpp"

namespace Attribute {

    class CmpOp {

        public:

            virtual bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) = 0;
            virtual ~CmpOp(){};
    };

    class IsEqCmpOp : public CmpOp {

        public:

            class IsEqVisitor : public boost::static_visitor<bool> {

                public:

                    template<typename T, typename U>
                    bool operator()(const T &lhs, const U &rhs) const {

                        assert(false);
                    }

                    template<typename T>
                    bool operator()(const T &lhs, const T &rhs) const {

                        return (lhs == rhs);
                    }

            };

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return boost::apply_visitor(IsEqVisitor(), lhs_data, rhs_data);
            };
    };

    class IsNotEqCmpOp : public CmpOp {

        public:

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return !boost::apply_visitor(IsEqCmpOp::IsEqVisitor(), lhs_data, rhs_data);
            };
    };

    class IsLessCmpOp : public CmpOp {

        public:

            class IsLessVisitor : public boost::static_visitor<bool> {

                public:

                    template<typename T, typename U>
                    bool operator()(const T &lhs, const U &rhs) const {

                        assert(false);
                    }

                    template<typename T>
                    bool operator()(const T &lhs, const T &rhs) const {

                        return (lhs < rhs);
                    }

            };

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return boost::apply_visitor(IsLessVisitor(), lhs_data, rhs_data);
            };
    };

    class IsLessOrEqCmpOp : public CmpOp {

        public:

            class IsLessOrEqVisitor : public boost::static_visitor<bool> {

                public:

                    template<typename T, typename U>
                    bool operator()(const T &lhs, const U &rhs) const {

                        assert(false);
                    }

                    template<typename T>
                    bool operator()(const T &lhs, const T &rhs) const {

                        return (lhs <= rhs);
                    }

            };

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return boost::apply_visitor(IsLessOrEqVisitor(), lhs_data, rhs_data);
            };
    };

    class IsGreaterCmpOp : public CmpOp {

        public:

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return !boost::apply_visitor(IsLessOrEqCmpOp::IsLessOrEqVisitor(), lhs_data, rhs_data);
            };
    };

    class IsGreaterOrEqCmpOp : public CmpOp {

        public:

            bool operator()(const Attribute::Value & lhs, const Attribute::Value & rhs) {

                auto lhs_data = lhs.data();
                auto rhs_data = rhs.data();

                return !boost::apply_visitor(IsLessCmpOp::IsLessVisitor(), lhs_data, rhs_data);
            };
    };
}

#endif
