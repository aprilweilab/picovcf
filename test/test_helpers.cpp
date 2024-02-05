#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getVCF42_EXAMPLE_FILE();

constexpr size_t DEFAULT_BUFFER_SIZE = 256;

/**
 * Utility class for reading/writing serialized data to memory instead
 * of disk. Just makes testing easier.
 */
class InMemBuffer : public std::streambuf
{
public:
    explicit InMemBuffer(size_t size)
            : m_buffer(new char[size]) {
        reset(size);
    }

    void reset(size_t size) {
        m_bufSize = size;
        memset(&m_buffer[0], 0, size);
        setp(m_buffer, m_buffer + m_bufSize);
        setg(m_buffer, m_buffer, m_buffer + m_bufSize);
    }

    size_t bytesWritten() const {
        return static_cast<size_t>(pptr() - pbase());
    }

    ~InMemBuffer() override {
        delete [] m_buffer;
    }

    char *m_buffer;
    size_t m_bufSize;
};

TEST(Helpers, Split) {
    std::string x = "a,b,c,defg";
    std::vector<std::string> tokens;
    picovcf::picovcf_split(x, ',', tokens);
    ASSERT_EQ(tokens, std::vector<std::string>({"a", "b", "c", "defg"}));

    std::string y = ",,,";
    tokens.clear();
    picovcf::picovcf_split(y, ',', tokens);
    ASSERT_EQ(tokens, std::vector<std::string>({"", "", "", ""}));
}


TEST(Helper, SerializeAlleleBits) {
    using SampleList = std::vector<IndexT>;
    InMemBuffer buffer(DEFAULT_BUFFER_SIZE);
    std::ostream outStream(&buffer);

    writeAllelesAsOnes(outStream, 1, {1, 1, 1, 1, 1, 1, 1, 1});
    ASSERT_EQ((uint8_t)buffer.m_buffer[0], 0xFF);
    ASSERT_EQ((uint8_t)buffer.m_buffer[1], 0x00);
    auto sampleSet = getSamplesWithAlt((const uint8_t*)&buffer.m_buffer[0], 8);
    ASSERT_EQ(sampleSet, SampleList({0, 1, 2, 3, 4, 5, 6, 7}));

    buffer.reset(DEFAULT_BUFFER_SIZE);
    writeAllelesAsOnes(outStream, 2, {2, 2, 2, 0, 1, 2, 2, 2, 2, 0});
    ASSERT_EQ((uint8_t)buffer.m_buffer[0], 0xE7);
    ASSERT_EQ((uint8_t)buffer.m_buffer[1], 0x80);
    sampleSet = getSamplesWithAlt((const uint8_t*)&buffer.m_buffer[0], 10);
    ASSERT_EQ(sampleSet, SampleList({0, 1, 2, 5, 6, 7, 8}));

    buffer.reset(DEFAULT_BUFFER_SIZE);
    writeAllelesAsOnes(outStream, 3, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3});
    ASSERT_EQ((uint8_t)buffer.m_buffer[0], 0x00);
    ASSERT_EQ((uint8_t)buffer.m_buffer[1], 0x01);
    sampleSet = getSamplesWithAlt((const uint8_t*)&buffer.m_buffer[0], 16);
    ASSERT_EQ(sampleSet, SampleList({15}));
}

TEST(Helper, BufferedReader) {
    std::string lineBuffer;
    std::vector<std::string> baseline;
    std::ifstream regularStream(getVCF42_EXAMPLE_FILE());
    size_t firstLineLen = 0;
    while (std::getline(regularStream, lineBuffer)) {
        if (firstLineLen == 0) {
            firstLineLen = lineBuffer.size();
        }
        baseline.push_back(lineBuffer);
    }
    ASSERT_EQ(baseline.size(), 25);

    // We use the line length + 1 as the buffer size to regression test a bug where we
    // concatenated lines together if the newline was exactly on the buffer boundary.
    BufferedReader br(getVCF42_EXAMPLE_FILE(), firstLineLen+1);
    size_t lineNo = 0;
    while (br.readline(lineBuffer)) {
        ASSERT_EQ(baseline.at(lineNo), lineBuffer);
        lineNo++;
    }
    ASSERT_EQ(baseline.size(), 25);
}

#if VCF_GZ_SUPPORT

extern const std::string getMSPRIME_EXAMPLE_FILE();
extern const std::string getMSPRIMEGZ_EXAMPLE_FILE();

TEST(Helper, ZBufferedReader) {
    std::string lineBuffer;
    std::vector<std::string> baseline;
    std::ifstream regularStream(getMSPRIME_EXAMPLE_FILE());
    while (std::getline(regularStream, lineBuffer)) {
        baseline.push_back(lineBuffer);
    }
    ASSERT_EQ(baseline.size(), 10);

    ZBufferedReader zbr(getMSPRIMEGZ_EXAMPLE_FILE(), 128);
    size_t lineNo = 0;
    while (zbr.readline(lineBuffer) > 0) {
        ASSERT_EQ(baseline.at(lineNo), lineBuffer);
        lineNo++;
    }
    ASSERT_EQ(lineNo, 10);
}
#endif