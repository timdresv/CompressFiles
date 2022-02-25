#include <iostream>
#include <windows.h>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <ctime>

class Compression {
private:
    std::vector<std::string> algorithms = {"SHF", "HUFF", "RLE", "BWT"};
    std::string dir_name, in_file_name, out_file_name;
    std::string mode, algorithm;
    HANDLE in_file, out_file;
    uint64_t in_buf_size, out_buf_size;
    BYTE *in_buffer, *out_buffer;
    clock_t time_alg;

    void Open() {
        SetCurrentDirectory(dir_name.c_str());
        in_file = CreateFile(in_file_name.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                             FILE_ATTRIBUTE_NORMAL, nullptr);
        DWORD err = GetLastError();
        if (err != 0)
            throw "Can't open file";
        in_buf_size = GetFileSize(in_file, nullptr);
        HANDLE in_file_map = CreateFileMapping(in_file, nullptr, PAGE_READONLY, 0, in_buf_size, in_file_name.c_str());
        in_buffer = (BYTE *) MapViewOfFile(in_file_map, FILE_MAP_READ, 0, 0, in_buf_size);
        CloseHandle(in_file_map);
    }

    void Create() {
        out_file = CreateFile(out_file_name.c_str(), GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ, nullptr, CREATE_NEW,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
        DWORD err = GetLastError();
        if (err != 0)
            throw "Can't create file";
        HANDLE out_file_map = CreateFileMapping(out_file, nullptr, PAGE_READWRITE, 0, out_buf_size,
                                                out_file_name.c_str());
        out_buffer = (BYTE *) MapViewOfFile(out_file_map, FILE_MAP_ALL_ACCESS, 0, 0, out_buf_size);
        CloseHandle(out_file_map);
    }

    void Close() {
        UnmapViewOfFile(in_buffer);
        UnmapViewOfFile(out_buffer);
        CloseHandle(in_file);
        CloseHandle(out_file);
        DeleteFile(in_file_name.c_str());
    }

    void Fail() {
        UnmapViewOfFile(in_buffer);
        UnmapViewOfFile(out_buffer);
        CloseHandle(in_file);
        CloseHandle(out_file);
        DeleteFile(out_file_name.c_str());
        throw "Failed to compress file";
    }


    void CompressMap(uint64_t code_size, int map_size, std::unordered_map<char, std::string> map_codes) {
        out_buf_size = in_buf_size;
        Create();
        uint64_t c = 0;
        if (7 > out_buf_size)
            Fail();
        out_buffer[c++] = distance(algorithms.begin(),
                                   find(algorithms.begin(), algorithms.end(), algorithm));//algorithm code(1B)
        uint64_t in_buf_size_temp = in_buf_size;
        for (int i = 0; i < 4; i++) {
            out_buffer[c++] = (BYTE) (in_buf_size_temp / ((uint64_t) (1) << 8 * (3 - i)));//size(4B)
            in_buf_size_temp %= ((uint64_t) (1) << 8 * (3 - i));
        }
        out_buffer[c++] = (BYTE) (map_codes.size() / 256);//map size(2B)
        out_buffer[c++] = (BYTE) (map_codes.size() % 256);
        BYTE b, two;
        for (auto &ch : map_codes) {
            if (c + 2 > out_buf_size)
                Fail();
            out_buffer[c++] = ch.first;
            out_buffer[c++] = (BYTE) ch.second.length();
            b = 0;
            two = 1 << (ch.second.length() - 1) % 8;
            for (size_t i = 0; i < ch.second.length(); i++) {
                if (ch.second[i] == '1') b += two;
                two /= 2;
                if (two == 0) {
                    if (c + 1 > out_buf_size)
                        Fail();
                    out_buffer[c++] = b;
                    b = 0;
                    two = 128;
                }
            }
        }
        if (c + 1 > out_buf_size)
            Fail();
        out_buffer[c++] = (code_size - 1) % 8;
        b = 0;
        two = 1 << (code_size - 1) % 8;
        std::string code;
        for (uint64_t i = 0; i < in_buf_size; i++) {
            code = map_codes[in_buffer[i]];
            for (size_t j = 0; j < code.length(); j++) {
                if (code[j] == '1') b += two;
                two /= 2;
                if (two == 0) {
                    if (c + 1 > out_buf_size)
                        Fail();
                    out_buffer[c++] = b;
                    b = 0;
                    two = 128;
                }
            }
        }
        out_buf_size = c;
        UnmapViewOfFile(out_buffer);
        SetFilePointer(out_file, out_buf_size % ((uint64_t) (1) << 32), (PLONG) (out_buf_size / ((uint64_t) (1) << 32)),
                       FILE_BEGIN);
        SetEndOfFile(out_file);
    }

    void DecompressMap() {
        uint64_t c = 1;
        out_buf_size = 0;
        for (int i = 0; i < 4; i++) {
            out_buf_size += in_buffer[c++] * ((uint64_t) 1 << 8 * (3 - i));
        }
        Create();
        int map_length = in_buffer[c++] * 256;
        map_length += in_buffer[c++];
        char ch;
        BYTE b, two;
        BYTE code_length;
        std::string code;
        std::unordered_map<std::string, char> map_char;
        for (int i = 0; i < map_length; i++) {
            ch = in_buffer[c++];
            code = "";
            code_length = in_buffer[c++];
            two = 1 << (code_length - 1) % 8;
            for (uint64_t j = c; c < j + (code_length - 1) / 8 + 1; c++) {
                b = in_buffer[c];
                while (two > 0) {
                    if (b / two == 1) {
                        code += "1";
                    } else {
                        code += "0";
                    }
                    b %= two;
                    two /= 2;
                }
                two = 128;
            }
            map_char.insert({code, ch});
        }
        uint64_t c2 = 0;
        two = 1 << in_buffer[c++];
        code = "";
        while (c < in_buf_size) {
            b = in_buffer[c++];
            while (two > 0) {
                if (b / two == 1) {
                    code += "1";
                } else {
                    code += "0";
                }
                b %= two;
                two /= 2;
                if (map_char.contains(code)) {
                    out_buffer[c2++] = map_char[code];
                    code = "";
                }
            }
            two = 128;
        }
    }

    std::unordered_map<char, uint64_t> FrequencyAnalisis() {
        std::unordered_map<char, uint64_t> map_frequency;
        for (uint64_t i = 0; i < in_buf_size; i++) {
            if (map_frequency.contains(in_buffer[i])) {
                map_frequency[in_buffer[i]]++;
            } else {
                map_frequency.insert({in_buffer[i], 1});
            }
        }
        return map_frequency;
    }


    struct Character {
        char ch;
        uint64_t count;
        std::string code;

        Character(char ch, uint64_t count, std::string code) : ch(ch), count(count), code(code) {}

        bool operator<(const Character &cmp) const {
            return count < cmp.count;
        }
    };

    std::vector<Character> AlgShenFano(std::vector<Character> vect, uint64_t sum) {
        if (vect.size() <= 1) return vect;
        uint64_t s = 0;
        size_t slice = 0;
        for (size_t i = 0; i < vect.size(); i++) {
            s += vect[i].count;
            if (s >= sum / 2 and i != 0) {
                s -= vect[i].count;
                slice = i;
                break;
            }
            vect[i].code += "0";
        }
        for (size_t i = slice; i < vect.size(); i++) {
            vect[i].code += "1";
        }
        std::vector<Character> vres = AlgShenFano(std::vector<Character>(vect.begin(), vect.begin() + slice), s);
        std::vector<Character> vtemp = AlgShenFano(std::vector<Character>(vect.begin() + slice, vect.end()), sum - s);
        vres.insert(vres.end(), vtemp.begin(), vtemp.end());
        return vres;
    }

    void CompressSHF() {
        std::unordered_map<char, uint64_t> map_frequency = FrequencyAnalisis();
        std::vector<Character> vect;
        for (auto &key : map_frequency) {
            vect.push_back(Character(key.first, key.second, ""));
        }
        sort(vect.begin(), vect.end());
        vect = AlgShenFano(vect, in_buf_size);
        uint64_t code_size = 0;
        int map_size = 0;
        std::unordered_map<char, std::string> map_codes;
        for (size_t i = 0; i < vect.size(); i++) {
            map_codes.insert({vect[i].ch, vect[i].code});
            code_size += vect[i].code.length() * vect[i].count;
            map_size += (((int) vect[i].code.length() - 1) / 8 + 1) + 2;
        }
        CompressMap(code_size, map_size, map_codes);
    }


    struct Node {
        char ch;
        uint64_t count;
        std::string code;
        Node *left, *right;

        Node(char ch, uint64_t count, std::string code) : ch(ch), count(count), code(code), left(nullptr),
                                                          right(nullptr) {}

        bool operator<(const Node &cmp) const {
            return count < cmp.count;
        }
    };

    void TreeTraversal(Node *tree, char ch) {
        if (tree->right == nullptr || tree->left == nullptr) {
            tree->code += ch;
        } else {
            TreeTraversal(tree->right, ch);
            TreeTraversal(tree->left, ch);
        }
    }

    void TreeToMap(Node *tree, uint64_t &code_size, int &map_size, std::unordered_map<char, std::string> &map_codes) {
        if (tree->right == nullptr || tree->left == nullptr) {
            std::reverse(tree->code.begin(), tree->code.end());
            map_codes.insert({tree->ch, tree->code});
            code_size += tree->code.length() * tree->count;
            map_size += (((int) tree->code.length() - 1) / 8 + 1) + 2;
        } else {
            TreeToMap(tree->right, code_size, map_size, map_codes);
            TreeToMap(tree->left, code_size, map_size, map_codes);
        }
    }

    void CompressHUFF() {
        std::unordered_map<char, uint64_t> map_frequency = FrequencyAnalisis();
        std::list<Node> lst;
        for (auto &key : map_frequency) {
            lst.push_back(Node(key.first, key.second, ""));
        }
        lst.sort();
        std::list<Node> tree_list;
        bool inserted;
        while (lst.size() >= 2) {
            Node *left = &(*lst.begin());
            Node *right = &(*(++lst.begin()));
            Node *parent = new Node(' ', left->count + right->count, "");
            TreeTraversal(left, '0');
            TreeTraversal(right, '1');
            tree_list.push_front(*left);
            tree_list.push_front(*right);
            parent->left = &(*tree_list.begin());
            parent->right = &(*(++tree_list.begin()));
            lst.erase(lst.erase(lst.begin()));
            inserted = false;
            for (auto iter = lst.begin(); iter != lst.end(); iter++) {
                if (parent->count <= iter->count) {
                    lst.insert(iter, *parent);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                lst.push_back(*parent);
            }
        }
        uint64_t code_size = 0;
        int map_size = 0;
        std::unordered_map<char, std::string> map_codes;
        TreeToMap(&lst.front(), code_size, map_size, map_codes);
        CompressMap(code_size, map_size, map_codes);
    }


    void CompressRLE() {
        out_buf_size = in_buf_size;
        Create();
        uint64_t c = 0;
        if (5 > out_buf_size)
            Fail();
        out_buffer[c++] = distance(algorithms.begin(),
                                   find(algorithms.begin(), algorithms.end(), algorithm));//algorithm code(1B)
        uint64_t in_buf_size_temp = in_buf_size;
        for (int i = 0; i < 4; i++) {
            out_buffer[c++] = (BYTE) (in_buf_size_temp / ((uint64_t) (1) << 8 * (3 - i)));//size(4B)
            in_buf_size_temp %= ((uint64_t) (1) << 8 * (3 - i));
        }
        int counter = 1;
        short md = 0;
        for (uint64_t i = 1; i < in_buf_size; i++) {
            if (in_buffer[i] == in_buffer[i - 1]) {
                if (md == 1) {
                    counter++;
                } else if (md == -1) {
                    //write counter-1 chars from in_buffer[i-counter] to in_buffer[i-2]
                    if (c + counter > out_buf_size)
                        Fail();
                    out_buffer[c++] = 128 + counter - 1;
                    for (int j = i - counter; j <= i - 2; j++) {
                        out_buffer[c++] = in_buffer[j];
                    }
                    md = 1;
                    counter = 2;
                } else {
                    md = 1;
                    counter++;
                }
            } else {
                if (md == -1) {
                    counter++;
                } else if (md == 1) {
                    //write counter chars in_buffer[i-1]
                    if (c + 2 > out_buf_size)
                        Fail();
                    out_buffer[c++] = counter;
                    out_buffer[c++] = in_buffer[i - 1];
                    md = 0;
                    counter = 1;
                } else {
                    md = -1;
                    counter++;
                }
            }
            if (counter >= 127 or i == in_buf_size - 1) {
                if (md == 1) {
                    //write counter chars in_buffer[i-1]
                    if (c + 2 > out_buf_size)
                        Fail();
                    out_buffer[c++] = counter;
                    out_buffer[c++] = in_buffer[i - 1];
                    md = 0;
                    counter = 0;
                } else if (md == -1) {
                    //write counter-1 chars from in_buffer[i-counter] to in_buffer[i-2]
                    if (c + counter + 1 > out_buf_size)
                        Fail();
                    out_buffer[c++] = 128 + counter;
                    for (int j = i - counter + 1; j <= i; j++) {
                        out_buffer[c++] = in_buffer[j];
                    }
                    md = 0;
                    counter = 0;
                }
            }
        }
        out_buf_size = c;
        UnmapViewOfFile(out_buffer);
        SetFilePointer(out_file, out_buf_size % ((uint64_t) (1) << 32), (PLONG) (out_buf_size / ((uint64_t) (1) << 32)),
                       FILE_BEGIN);
        SetEndOfFile(out_file);
    }

    void DecompressRLE() {
        uint64_t c = 1, c2 = 0;
        out_buf_size = 0;
        for (int i = 0; i < 4; i++) {
            out_buf_size += in_buffer[c++] * ((uint64_t) 1 << 8 * (3 - i));
        }
        Create();
        while (c < in_buf_size && c2 < out_buf_size) {
            int n = in_buffer[c++];
            if (n < 128) {
                for (int i = 0; i < n; i++) {
                    out_buffer[c2++] = in_buffer[c];
                }
                c++;
            } else {
                for (int i = 0; i < n - 128; i++) {
                    out_buffer[c2++] = in_buffer[c++];
                }
            }
        }
    }


    void ModificationBWT() {
        int block_size = 256, out_index = 0;
        out_buf_size = in_buf_size + in_buf_size / block_size + 2;
        Create();
        uint64_t c = 0;
        out_buffer[c++] = distance(algorithms.begin(),
                                   find(algorithms.begin(), algorithms.end(), algorithm));//algorithm code(1B)
        std::vector<std::string> block(block_size);
        for (uint64_t i = 0; i < in_buf_size; i += block_size) {
            if (block_size > in_buf_size - i) {
                block_size = in_buf_size - i;
                block.resize(block_size);
            }
            for (int j = 0; j < block_size; j++) {
                for (int l = j; l < block_size; l++) {
                    block[j] += in_buffer[i + l];
                }
                for (int l = 0; l < j; l++) {
                    block[j] += in_buffer[i + l];
                }
            }
            block[0] += "0";
            sort(block.begin(), block.end());
            for (int j = 0; j < block_size; j++) {
                out_buffer[c++] = block[j][block_size - 1];
                if (block[j].length() != block_size) {
                    out_index = j;
                }
                block[j] = "";
            }
            out_buffer[c++] = out_index;
        }
    }

    void DemodificationBWT() {
        int block_size = 256, out_index, ind = 0;
        out_buf_size = in_buf_size;
        Create();
        uint64_t c = 0;
        std::string block, block_sort;
        std::vector<int> key(block_size);
        for (uint64_t i = 1; i < in_buf_size; i += block_size + 1) {
            block_size = (block_size < in_buf_size - 1 - i) ? block_size : in_buf_size - 1 - i;
            block = "";
            for (int j = 0; j < block_size; j++) {
                block += in_buffer[i + j];
            }
            out_index = in_buffer[i + block_size];
            block_sort = block;
            sort(block_sort.begin(), block_sort.end());
            ind = 0;
            for (int j = 0; j < block_size; j++) {
                if (block_sort[(j != 0) ? j - 1 : j] != block_sort[j])
                    ind = 0;
                while (block[ind] != block_sort[j]) {
                    if (ind == block_size - 1)
                        ind = 0;
                    else
                        ind += 1;
                }
                key[j] = ind;
                if (ind == block_size - 1)
                    ind = 0;
                else
                    ind += 1;
            }
            ind = out_index;
            do {
                ind = key[ind];
                out_buffer[c++] = block[ind];
            } while (ind != out_index);
        }
        out_buf_size = c;
        UnmapViewOfFile(out_buffer);
        SetFilePointer(out_file, out_buf_size % ((uint64_t) (1) << 32), (PLONG) (out_buf_size / ((uint64_t) (1) << 32)),
                       FILE_BEGIN);
        SetEndOfFile(out_file);
    }

public:
    Compression(std::string file_name, std::string md, std::string alg) {
        clock_t start_time;
        size_t ind = file_name.find_last_of('\\');
        dir_name = file_name.substr(0, ind);
        in_file_name = file_name.substr(ind + 1);
        mode = md;
        if (mode == "compress") {
            out_file_name = in_file_name + ".cmp";
            Open();
            algorithm = alg;
            start_time = clock();
            if (algorithm == "SHF") {
                CompressSHF();
            } else if (algorithm == "HUFF") {
                CompressHUFF();
            } else if (algorithm == "RLE") {
                CompressRLE();
            } else if (algorithm == "BWT") {
                ModificationBWT();
            } else {
                Close();
                throw "Wrong algorithm";
            }
        } else if (mode == "decompress") {
            ind = in_file_name.find_last_of('.');
            out_file_name = in_file_name.substr(0, ind);
            Open();
            algorithm = algorithms[in_buffer[0]];
            start_time = clock();
            if (algorithm == "SHF" || algorithm == "HUFF") {
                DecompressMap();
            } else if (algorithm == "RLE") {
                DecompressRLE();
            } else if (algorithm == "BWT") {
                DemodificationBWT();
            } else {
                Close();
                throw "Wrong file's extension";
            }
        } else {
            throw "Wrong mode";
        }
        clock_t end_time = clock();
        time_alg = end_time - start_time;
        Close();
    }

    void PrintInformation() {
        if (mode == "compress") {
            std::cout << std::endl;
            std::cout << "Algorithm " << algorithm << ": " << std::endl;
            std::cout << "\tSource file size: " << in_buf_size / 1024.0 << " Kb" << std::endl;
            std::cout << "\tCompressed file size: " << out_buf_size / 1024.0 << " Kb" << std::endl;
            std::cout << "\tTime compressing: " << time_alg / 1000.0 << " sec" << std::endl;
            std::cout << std::endl;
        } else if (mode == "decompress") {
            std::cout << std::endl;
            std::cout << "Algorithm " << algorithm << ": " << std::endl;
            std::cout << "\tCompressed file size: " << in_buf_size / 1024.0 << " Kb" << std::endl;
            std::cout << "\tDecompressed file size: " << out_buf_size / 1024.0 << " Kb" << std::endl;
            std::cout << "\tTime decompressing: " << time_alg / 1000.0 << " sec" << std::endl;
            std::cout << std::endl;
        }
    }
};

std::string OpenFileDialog() {
    OPENFILENAME ofn;
    TCHAR szFile[260] = {0};
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = nullptr;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = nullptr;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn) == TRUE) {
        return ofn.lpstrFile;
    } else {
        return nullptr;
    }
}

int main() {
    std::string input;
    std::string mode, alg;
    while (true) {
        std::cout << "Enter encoding mode (compress/decompress): ";
        std::cin >> mode;
        if (mode == "compress") {
            std::cout << "Enter compression algorithm (SHF/HUFF/RLE/BWT): ";
            std::cin >> alg;
        }
        input = OpenFileDialog();
        try {
            Compression cmp(input, mode, alg);
            cmp.PrintInformation();
        }
        catch (const char *msg) {
            std::cerr << msg << std::endl;
        }
    }
    return 0;
}