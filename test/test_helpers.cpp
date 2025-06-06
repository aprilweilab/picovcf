#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <limits>
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
    std::vector<uint8_t> buffer;

    writeAllelesAsOnes(buffer, {0, 1, 2, 3, 4, 5, 6, 7}, 8);
    ASSERT_EQ((uint8_t)buffer.at(0), 0xFF);
    auto sampleSet = getSamplesWithAlt((const uint8_t*)buffer.data(), 8);
    ASSERT_EQ(sampleSet, IGDSampleList({0, 1, 2, 3, 4, 5, 6, 7}));

	buffer.clear();
    writeAllelesAsOnes(buffer, {0, 1, 2, 5, 6, 7, 8}, 10);
    ASSERT_EQ((uint8_t)buffer.at(0), 0xE7);
    ASSERT_EQ((uint8_t)buffer.at(1), 0x80);
    sampleSet = getSamplesWithAlt((const uint8_t*)buffer.data(), 10);
    ASSERT_EQ(sampleSet, IGDSampleList({0, 1, 2, 5, 6, 7, 8}));

    buffer.clear();
    writeAllelesAsOnes(buffer, {15}, 16);
    ASSERT_EQ((uint8_t)buffer.at(0), 0x00);
    ASSERT_EQ((uint8_t)buffer.at(1), 0x01);
    sampleSet = getSamplesWithAlt((const uint8_t*)buffer.data(), 16);
    ASSERT_EQ(sampleSet, IGDSampleList({15}));
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

TEST(Helper, IGDAllele) {
    auto a1 = IGDAllele("", 0);
    ASSERT_FALSE(a1.isLong());
    ASSERT_EQ(a1.getString(), "");

    auto a2 = IGDAllele("A", 0);
    ASSERT_FALSE(a2.isLong());
    ASSERT_EQ(a2.getString(), "A");

    auto a3 = IGDAllele("ACG", 0);
    ASSERT_FALSE(a3.isLong());
    ASSERT_EQ(a3.getString(), "ACG");

    auto a4 = IGDAllele("ACGT", 0);
    ASSERT_FALSE(a4.isLong());
    ASSERT_EQ(a4.getString(), "ACGT");

    auto a5 = IGDAllele("ACGTC", 0);
    ASSERT_TRUE(a5.isLong());
    ASSERT_EQ(a5.getLongIndex(), 0);

    EXPECT_THROW(IGDAllele("ACGTC", ((size_t)std::numeric_limits<uint32_t>::max()) + 500), MalformedFile);

    std::string fake = "ACGT";
    fake[3] = 0x80; // Bad string value
    EXPECT_THROW(IGDAllele(fake, 0), MalformedFile);

    size_t largestIndex = std::numeric_limits<uint32_t>::max() / 2;
    auto a6 = IGDAllele("ACGTC", largestIndex);
    ASSERT_TRUE(a6.isLong());
    ASSERT_EQ(a6.getLongIndex(), largestIndex);

}

TEST(Helper, StructuredMetadata) {
    std::map<std::string, std::string> result;
    std::map<std::string, std::string> expected;
    // Doesn't start/end with <...>
    EXPECT_THROW(
        result = picovcf_parse_structured_meta("no <> characters"),
        MalformedFile);
    // Unterminated quote
    EXPECT_THROW(
        result = picovcf_parse_structured_meta("<A=\"blah>"),
        MalformedFile);

    result = picovcf_parse_structured_meta("<>");
    expected.clear();
    EXPECT_EQ(result, expected);

    result = picovcf_parse_structured_meta("<A=B,C=D>");
    expected = {{"A", "B"}, {"C", "D"}};
    EXPECT_EQ(result, expected);
    result = picovcf_parse_structured_meta("<A=\"B\",C=\"D\">");
    EXPECT_EQ(result, expected);

    result = picovcf_parse_structured_meta("<A=\"this has = and , and > characters in a string\">");
    expected = {{"A", "this has = and , and > characters in a string"}};
    EXPECT_EQ(result, expected);
}