template<class T, int size1, int size2>
class vector2{

    private:
        T data[size1][size2];

        friend class proxy;

    class proxy { 
        vector2 &m_;
        int index1_;
        public:
            proxy(vector2 &m, int i1) : 
                m_(m), index1_(i1) {}
            T &operator[](int index2) { 
                return m_.data[index1_][index2];
            }
    };
    
    public:
        proxy operator[](int index) {
            return proxy(*this, index);
        }
};

template<class T, int size1, int size2, int size3>
class vector3{

    private:
        T data[size1][size1][size3];

        friend class proxy;
        friend class proxy2;

    class proxy { 
        vector3 &m_;
        int index1_, index2_;
        public:
            proxy(vector3 &m, int i1, int i2) : 
                m_(m), index1_(i1), index2_(i2) {}
            T &operator[](int index3) { 
                return m_.data[index1_][index2_][index3];
            }
    };
    
    class proxy2 { 
        vector3 &m_;
        int index_;
        public:
            proxy2(vector3 &m, int d) : m_(m), index_(d) { }
            proxy operator[](int index2) { 
                return proxy(m_, index_, index2);
            }
    };
    
    public:
        proxy2 operator[](int index) {
            return proxy2(*this, index);
        }
};
