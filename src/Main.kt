import java.awt.image.BufferedImage
import java.io.DataOutputStream
import java.io.File
import java.io.FileOutputStream
import java.io.IOException
import javax.imageio.ImageIO
import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sqrt


enum class DownsamplingMethod {
    REMOVE_LINES, AVERAGE, CLOSEST
}

fun getQuantizationMatrix(Q: Int): Array<IntArray> {
    var q = Q
    val standardQ = arrayOf(
        intArrayOf(16, 11, 10, 16, 24, 40, 51, 61),
        intArrayOf(12, 12, 14, 19, 26, 58, 60, 55),
        intArrayOf(14, 13, 16, 24, 40, 57, 69, 56),
        intArrayOf(14, 17, 22, 29, 51, 87, 80, 62),
        intArrayOf(18, 22, 37, 56, 68, 109, 103, 77),
        intArrayOf(24, 35, 55, 64, 81, 104, 113, 92),
        intArrayOf(49, 64, 78, 87, 103, 121, 120, 101),
        intArrayOf(72, 92, 95, 98, 112, 100, 103, 99)
    )
    var scale = 1.0
    if (q == 0) {
        q = 1
    }
    if (q < 50) {
        scale = 5000.0 / q
    } else {
        scale = (200.0 - 2.0 * q)
    }

    val quantizationMatrix = Array(8) { IntArray(8) }
    for (i in 0..7) {
        for (j in 0..7) {
            quantizationMatrix[i][j] = ((standardQ[i][j] * scale) + 50).toInt() / 100
            if (quantizationMatrix[i][j] == 0) {
                quantizationMatrix[i][j] = 1
            }
        }
    }

    //val quantizationMatrix = Array(8) { IntArray(8) }
    //if (Q < 50) {
    //    for (i in 0 until 8) {
    //        for (j in 0 until 8) {
    //            quantizationMatrix[i][j] = ((50 + Q) * standardQ[i][j] / 50).coerceAtLeast(1)
    //        }
    //    }
    //} else {
    //    for (i in 0 until 8) {
    //        for (j in 0 until 8) {
    //            quantizationMatrix[i][j] = (2 * (100 - Q) * standardQ[i][j] / 100).coerceAtLeast(1)
    //        }
    //    }
    //}

    return quantizationMatrix
}

fun rgbToYCbCr(rgb: Array<Array<IntArray>>): Array<Array<IntArray>> {
    val yCbCr = Array(rgb.size) { Array(rgb[0].size) { IntArray(3) } }
    for (i in rgb.indices) {
        for (j in rgb[0].indices) {
            val r = rgb[i][j][0].toDouble()
            val g = rgb[i][j][1].toDouble()
            val b = rgb[i][j][2].toDouble()

            val y = 0.299 * r + 0.587 * g + 0.114 * b
            val cb = 128 - 0.168736 * r - 0.331264 * g + 0.5 * b
            val cr = 128 + 0.5 * r - 0.418688 * g - 0.081312 * b

            yCbCr[i][j][0] = y.toInt().coerceIn(0, 255)
            yCbCr[i][j][1] = cb.toInt().coerceIn(0, 255)
            yCbCr[i][j][2] = cr.toInt().coerceIn(0, 255)
        }
    }
    return yCbCr
}

fun yCbCrToRgb(yCbCr: Array<Array<IntArray>>): Array<Array<IntArray>> {
    val rgb = Array(yCbCr.size) { Array(yCbCr[0].size) { IntArray(3) } }
    for (i in yCbCr.indices) {
        for (j in yCbCr[0].indices) {
            val y = yCbCr[i][j][0].toDouble()
            val cb = yCbCr[i][j][1].toDouble() - 128
            val cr = yCbCr[i][j][2].toDouble() - 128

            val r = (y + 1.402 * cr).toInt().coerceIn(0, 255)
            val g = (y - 0.344136 * cb - 0.714136 * cr).toInt().coerceIn(0, 255)
            val b = (y + 1.772 * cb).toInt().coerceIn(0, 255)

            rgb[i][j][0] = r
            rgb[i][j][1] = g
            rgb[i][j][2] = b
        }
    }
    return rgb
}

fun zigzag(matrix: Array<IntArray>): IntArray {
    val n = matrix.size
    val m = matrix[0].size
    val result = IntArray(n * m)
    var i = 0
    var j = 0
    for (k in result.indices) {
        result[k] = matrix[i][j]
        if ((i + j) % 2 == 0) { // moving up
            if (j == m - 1) i++
            else if (i == 0) j++
            else {
                i--; j++
            }
        } else { // moving down
            if (i == n - 1) j++
            else if (j == 0) i++
            else {
                i++; j--
            }
        }
    }
    return result
}

fun dezigzag(array: IntArray, rows: Int, cols: Int): Array<IntArray> {
    val result = Array(rows) { IntArray(cols) }
    var i = 0
    var j = 0
    for (k in array.indices) {
        result[i][j] = array[k]
        if ((i + j) % 2 == 0) { // moving up
            if (j == cols - 1) i++
            else if (i == 0) j++
            else {
                i--; j++
            }
        } else { // moving down
            if (i == rows - 1) j++
            else if (j == 0) i++
            else {
                i++; j--
            }
        }
    }
    return result
}

data class RLEElement(val value: Int, val count: Int)

fun rleEncode(data: IntArray): List<RLEElement> {
    if (data.isEmpty()) return emptyList()

    val encoded = mutableListOf<RLEElement>()
    var currentValue = data[0]
    var count = 1

    for (i in 1 until data.size) {
        if (data[i] == currentValue) {
            count++
        } else {
            encoded.add(RLEElement(currentValue, count))
            currentValue = data[i]
            count = 1
        }
    }
    encoded.add(RLEElement(currentValue, count)) // добавляем последнюю пару
    return encoded
}

fun rleDecode(encoded: List<RLEElement>): IntArray {
    val decoded = mutableListOf<Int>()
    for (element in encoded) {
        repeat(element.count) {
            decoded.add(element.value)
        }
    }
    return decoded.toIntArray()
}

fun encode3DArray(data: Array<Array<Array<IntArray>>>): List<List<List<List<RLEElement>>>> {
    return data.map { twoDArray ->
        twoDArray.map { oneDArray ->
            oneDArray.map { intArray ->
                rleEncode(intArray)
            }
        }
    }
}

fun decode3DArray(encoded: List<List<List<List<RLEElement>>>>): Array<Array<Array<IntArray>>> {
    return Array(encoded.size) { i ->
        Array(encoded[i].size) { j ->
            Array(encoded[i][j].size) { k ->
                rleDecode(encoded[i][j][k])
            }
        }
    }
}

fun downsample(matrix: Array<IntArray>, cx: Int, cy: Int, method: DownsamplingMethod): Array<IntArray> {
    val newWidth = matrix.size / cx
    val newHeight = matrix[0].size / cy
    val result = Array(newWidth) { IntArray(newHeight) }

    when (method) {
        DownsamplingMethod.REMOVE_LINES -> {
            for (i in 0 until newWidth) {
                for (j in 0 until newHeight) {
                    result[i][j] = matrix[i * cx][j * cy]
                }
            }
        }

        DownsamplingMethod.AVERAGE -> {
            for (i in 0 until newWidth) {
                for (j in 0 until newHeight) {
                    var sum = 0
                    for (x in 0 until cx) {
                        for (y in 0 until cy) {
                            sum += matrix[i * cx + x][j * cy + y]
                        }
                    }
                    result[i][j] = sum / (cx * cy)
                }
            }
        }

        DownsamplingMethod.CLOSEST -> {
            for (i in 0 until newWidth) {
                for (j in 0 until newHeight) {
                    var sum = 0
                    for (x in 0 until cx) {
                        for (y in 0 until cy) {
                            sum += matrix[i * cx + x][j * cy + y]
                        }
                    }
                    val avg = sum / (cx * cy)
                    var closestValue = matrix[i * cx][j * cy]
                    var minDiff = Int.MAX_VALUE
                    for (x in 0 until cx) {
                        for (y in 0 until cy) {
                            val diff = abs(matrix[i * cx + x][j * cy + y] - avg)
                            if (diff < minDiff) {
                                minDiff = diff
                                closestValue = matrix[i * cx + x][j * cy + y]
                            }
                        }
                    }
                    result[i][j] = closestValue
                }
            }
        }
    }
    return result
}

fun upsample(matrix: Array<IntArray>, cx: Int, cy: Int): Array<IntArray> {
    val newWidth = matrix.size * cx
    val newHeight = matrix[0].size * cy
    val result = Array(newWidth) { IntArray(newHeight) }

    for (i in matrix.indices) {
        for (j in matrix[0].indices) {
            for (x in 0 until cx) {
                for (y in 0 until cy) {
                    result[i * cx + x][j * cy + y] = matrix[i][j]
                }
            }
        }
    }
    return result
}

fun dct(matrix: Array<IntArray>): Array<IntArray> {
    val N = matrix.size
    val dctMatrix = Array(N) { IntArray(N) }
    val C = Array(N) { DoubleArray(N) }

    for (i in 0 until N) {
        for (j in 0 until N) {
            C[i][j] = cos((2 * i + 1) * j * PI / (2.0 * N))
        }
    }

    for (u in 0 until N) {
        for (v in 0 until N) {
            var sum = 0.0
            for (x in 0 until N) {
                for (y in 0 until N) {
                    sum += matrix[x][y] * C[x][u] * C[y][v]
                }
            }
            val cu = if (u == 0) 1 / sqrt(2.0) else 1.0
            val cv = if (v == 0) 1 / sqrt(2.0) else 1.0
            dctMatrix[u][v] = (0.25 * cu * cv * sum).toInt()
        }
    }
    return dctMatrix
}

fun idct(matrix: Array<IntArray>): Array<IntArray> {
    val N = matrix.size
    val idctMatrix = Array(N) { IntArray(N) }
    val C = Array(N) { DoubleArray(N) }

    for (i in 0 until N) {
        for (j in 0 until N) {
            C[i][j] = cos((2 * i + 1) * j * PI / (2.0 * N))
        }
    }

    for (x in 0 until N) {
        for (y in 0 until N) {
            var sum = 0.0
            for (u in 0 until N) {
                for (v in 0 until N) {
                    val cu = if (u == 0) 1 / sqrt(2.0) else 1.0
                    val cv = if (v == 0) 1 / sqrt(2.0) else 1.0
                    sum += cu * cv * matrix[u][v] * C[x][u] * C[y][v]
                }
            }
            idctMatrix[x][y] = (0.25 * sum).toInt()
        }
    }
    return idctMatrix
}

fun quantize(matrix: Array<IntArray>, quantizationTable: Array<IntArray>): Array<IntArray> {
    val result = Array(matrix.size) { IntArray(matrix[0].size) }
    for (i in matrix.indices) {
        for (j in matrix[0].indices) {
            result[i][j] = (matrix[i][j] / quantizationTable[i][j])
        }
    }
    return result
}

fun dequantize(matrix: Array<IntArray>, quantizationTable: Array<IntArray>): Array<IntArray> {
    val result = Array(matrix.size) { IntArray(matrix[0].size) }
    for (i in matrix.indices) {
        for (j in matrix[0].indices) {
            result[i][j] = matrix[i][j] * quantizationTable[i][j]
        }
    }
    return result
}

fun saveToFile(data: String, fileName: String) {
    val file = File(fileName)
    file.writeText(data)
}

fun write3DArrayToFile(array3D: Array<Array<IntArray>>, filename: String) {
    try {
        DataOutputStream(FileOutputStream(filename)).use { dos ->
            // Записываем размеры массива
//            dos.writeInt(array3D.length);
//            dos.writeInt(array3D[0].length);
//            dos.writeInt(array3D[0][0].length);

            // Записываем элементы массива
            for (matrix in array3D) {
                for (row in matrix) {
                    for (value in row) {
                        dos.writeByte(value)
                    }
                }
            }
        }
    } catch (e: IOException) {
        e.printStackTrace()
    }
}

fun write2DArrayToFile(array2D: Array<IntArray>, filename: String) {
    try {
        DataOutputStream(FileOutputStream(filename)).use { dos ->
            // Записываем размеры массива
            dos.writeInt(array2D.size)
            dos.writeInt(array2D[0].size)

            // Записываем элементы массива
            for (row in array2D) {
                for (value in row) {
                    dos.writeByte(value)
                }
            }
        }
    } catch (e: IOException) {
        e.printStackTrace()
    }
}

fun calculateCompressedSize(rleData: List<RLEElement>, bitsPerValue: Int, bitsPerCount: Int): Int {
    val bitsPerElement = bitsPerValue + bitsPerCount
    return rleData.size * bitsPerElement
}

fun rleEncodeArray(input: IntArray): IntArray {
    if (input.isEmpty()) return IntArray(0)

    val encodedList = mutableListOf<Int>()
    var currentNumber = input[0]
    var count = 1

    for (i in 1 until input.size) {
        if (input[i] == currentNumber) {
            count++
        } else {
            encodedList.add(currentNumber)
            encodedList.add(count)
            currentNumber = input[i]
            count = 1
        }
    }
    // Добавляем последнюю группу
    encodedList.add(currentNumber)
    encodedList.add(count)

    return encodedList.toIntArray()
}

fun rleDecodeArray(encoded: IntArray): IntArray {
    if (encoded.isEmpty()) return IntArray(0)

    val decodedList = mutableListOf<Int>()

    for (i in encoded.indices step 2) {
        val value = encoded[i]
        val count = encoded[i + 1]
        repeat(count) {
            decodedList.add(value)
        }
    }

    return decodedList.toIntArray()
}


fun main() {
    var p = 50
    //while (p <= 100) {
        val imageName = "万里长城" //万里长城 random purple
        // Загружаем изображение
        val image = ImageIO.read(File("D:\\ЛЭТИ\\АиСД\\Второй семестр\\Изображения\\$imageName.jpg"))
        val width = image.width
        val height = image.height
        val rgb = Array(width) { Array(height) { IntArray(3) } }
        for (x in 0 until width) {
            for (y in 0 until height) {
                val color = image.getRGB(x, y)
                rgb[x][y][0] = (color shr 16) and 0xFF
                rgb[x][y][1] = (color shr 8) and 0xFF
                rgb[x][y][2] = color and 0xFF
            }
        }
        write3DArrayToFile(rgb, "1.rgbFile.txt")

        // Преобразуем RGB в YCbCr
        val yCbCr = rgbToYCbCr(rgb)
        //val file1 = "D:\\ЛЭТИ\\АиСД\\Второй семестр\\Размеры\\yCbCrFile.txt"
        //val file2 = "D:\\ЛЭТИ\\АиСД\\Второй семестр\\Размеры\\downsampleFile.txt"
        //val file3 = "D:\\ЛЭТИ\\АиСД\\Второй семестр\\Размеры\\yCbCrFile.txt"
        //val file4 = "D:\\ЛЭТИ\\АиСД\\Второй семестр\\Размеры\\yCbCrFile.txt"
        write3DArrayToFile(yCbCr, "2.yCbCrFile.txt")
        //println("Size after YCbCr (bits): " + file1.length * 8)
        //val fileName1 = "1.YCbCr_$imageName.txt"
        //saveToFile(yCbCr.map { it.map { it[0] }.toIntArray() }.toTypedArray().toString() +
        //        yCbCr.map { it.map { it[1] }.toIntArray() }.toTypedArray().toString() +
        //        yCbCr.map { it.map { it[2] }.toIntArray() }.toTypedArray().toString(), fileName1)

        // Downsample Cb and Cr channels
        //val downsampledY = downsample(yCbCr.map { it.map { it[0] }.toIntArray() }.toTypedArray(), 2, 2, DownsamplingMethod.AVERAGE)
        val downsampledCb =
            downsample(yCbCr.map { it.map { it[1] }.toIntArray() }.toTypedArray(), 2, 2, DownsamplingMethod.AVERAGE)
        val downsampledCr =
            downsample(yCbCr.map { it.map { it[2] }.toIntArray() }.toTypedArray(), 2, 2, DownsamplingMethod.AVERAGE)

        write2DArrayToFile(yCbCr.map { it.map { it[0] }.toIntArray() }.toTypedArray(), "3.YChanel.txt")
        write2DArrayToFile(downsampledCb, "3.downsampledCb.txt")
        write2DArrayToFile(downsampledCr, "3.downsampledCr.txt")

        // Применяем DCT и квантование
        val blockSize = 8
        val quality = p // Уровень качества от 1 до 100
        val quantizationTable = getQuantizationMatrix(quality)
        var dctSize: String = ""
        var quantSize: String = ""
        val dctCoefficients =
            Array(3) { Array(width / blockSize) { Array(height / blockSize) { Array(blockSize) { IntArray(blockSize) } } } }
        val dctSizeCounterY = Array(width / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        val quantSizeCounterY = Array(width / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        var rleSizeCounterY = 0
        var rleSizeCounterCb = 0
        var rleSizeCounterCr = 0
        //DCT & quant только канала Y
        for (i in 0..<width / blockSize) {
            for (j in 0..<height / blockSize) {
                val block = Array(blockSize) { IntArray(blockSize) }
                for (x in 0..<blockSize) {
                    for (y in 0..<blockSize) {
                        block[x][y] = yCbCr[i * blockSize + x][j * blockSize + y][0]
                    }
                }

                val dctBlock = dct(block)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        dctSizeCounterY[i][x][y] = dctBlock[x][y]
                    }
                }
                val quantizedBlock = quantize(dctBlock, quantizationTable)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        quantSizeCounterY[i][x][y] = quantizedBlock[x][y]
                    }
                }
                val zigzagBlock = zigzag(quantizedBlock)
                //val rleBlock: List<RLEElement> = rleEncode(zigzagBlock) //переписать RLE
                val rleBlock = rleEncodeArray(zigzagBlock)
                rleSizeCounterY += rleBlock.size * 4
                // Store RLE encoded block
                dctCoefficients[0][i][j] = dezigzag(rleDecodeArray(rleBlock), blockSize, blockSize)
            }
        }

        write3DArrayToFile(dctSizeCounterY, "4.dctSizeY.txt")
        write3DArrayToFile(quantSizeCounterY, "5.quantSizeY.txt")

        val dctSizeCounterCb = Array(width / 2 / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        val quantSizeCounterCb = Array(width / 2 / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        //DCT & quant только канала Cb
        for (i in 0..<width / 2 / blockSize) {
            for (j in 0..<height / 2 / blockSize) {
                val block = Array(blockSize) { IntArray(blockSize) }
                for (x in 0..<blockSize) {
                    for (y in 0..<blockSize) {
                        block[x][y] = downsampledCb[i * blockSize + x][j * blockSize + y]
                    }
                }

                val dctBlock = dct(block)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        dctSizeCounterCb[i][x][y] = dctBlock[x][y]
                    }
                }
                val quantizedBlock = quantize(dctBlock, quantizationTable)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        quantSizeCounterCb[i][x][y] = quantizedBlock[x][y]
                    }
                }
                val zigzagBlock = zigzag(quantizedBlock)
                //val rleBlock: List<RLEElement> = rleEncode(zigzagBlock) //переписать RLE
                val rleBlock = rleEncodeArray(zigzagBlock)
                //println(dezigzag(rleDecode(rleBlock), blockSize, blockSize).toString())
                rleSizeCounterCb += rleBlock.size * 4

                // Store RLE encoded block
                dctCoefficients[1][i][j] = dezigzag(rleDecodeArray(rleBlock), blockSize, blockSize)
            }
        }

        write3DArrayToFile(dctSizeCounterCb, "4.dctSizeCb.txt")
        write3DArrayToFile(quantSizeCounterCb, "5.quantSizeCb.txt")

        val dctSizeCounterCr = Array(width / 2 / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        val quantSizeCounterCr = Array(width / 2 / blockSize) { Array(blockSize) { IntArray(blockSize) } }
        //DCT & quant только канала Cr
        for (i in 0..<width / 2 / blockSize) {
            for (j in 0..<height / 2 / blockSize) {
                val block = Array(blockSize) { IntArray(blockSize) }
                for (x in 0..<blockSize) {
                    for (y in 0..<blockSize) {
                        block[x][y] = downsampledCr[i * blockSize + x][j * blockSize + y]
                    }
                }

                val dctBlock = dct(block)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        dctSizeCounterCr[i][x][y] = dctBlock[x][y]
                    }
                }
                val quantizedBlock = quantize(dctBlock, quantizationTable)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        quantSizeCounterCr[i][x][y] = quantizedBlock[x][y]
                    }
                }
                val zigzagBlock = zigzag(quantizedBlock)
                //val rleBlock: List<RLEElement> = rleEncode(zigzagBlock) //переписать RLE
                val rleBlock = rleEncodeArray(zigzagBlock)
                //println(dezigzag(rleDecode(rleBlock), blockSize, blockSize).toString())
                rleSizeCounterCr += rleBlock.size * 4

                // Store RLE encoded block
                dctCoefficients[2][i][j] = dezigzag(rleDecodeArray(rleBlock), blockSize, blockSize)
            }
        }

        println("${rleSizeCounterY + rleSizeCounterCb + rleSizeCounterCr}")
        write3DArrayToFile(dctSizeCounterCr, "4.dctSizeCr.txt")
        write3DArrayToFile(quantSizeCounterCr, "5.quantSizeCr.txt")


        //val fileName3 = "3.dct_$imageName.txt"
        //saveToFile(dctSize, fileName3)
//
        //val fileName4 = "4.quant_$imageName.txt"
        //saveToFile(quantSize, fileName4)

        //val encodedRLE = encode3DArray(dctCoefficients)
        //val decodeRLE = decode3DArray(encodedRLE)


        // Применяем деквантование и обратное DCT к каналу Y
        val yCbCrDecoded = Array(width) { Array(height) { IntArray(3) } }
        for (i in 0 until width / blockSize) {
            for (j in 0 until height / blockSize) {
                val quantizedBlock = dctCoefficients[0][i][j]
                val dequantizedBlock = dequantize(quantizedBlock, quantizationTable)
                val idctBlock = idct(dequantizedBlock)
                for (x in 0 until blockSize) {
                    for (y in 0 until blockSize) {
                        yCbCrDecoded[i * blockSize + x][j * blockSize + y][0] = idctBlock[x][y]
                    }
                }
            }
        }

        // Применяем деквантование и обратное DCT к каналу Cb
        //var toUpsacleCb: Array<IntArray>
        for (i in 0 until width / 2 / blockSize) {
            for (j in 0 until height / 2 / blockSize) {
                val quantizedBlock = dctCoefficients[1][i][j]
                val dequantizedBlock = dequantize(quantizedBlock, quantizationTable)
                val idctBlock = upsample(idct(dequantizedBlock), 2, 2)
                for (x in 0 until blockSize * 2) {
                    for (y in 0 until blockSize * 2) {
                        yCbCrDecoded[i * blockSize + x][j * blockSize + y][1] = idctBlock[x][y]
                    }
                }
            }
        }

        // Применяем деквантование и обратное DCT к каналу Cr
        for (i in 0 until width / 2 / blockSize) {
            for (j in 0 until height / 2 / blockSize) {
                val quantizedBlock = dctCoefficients[2][i][j]
                val dequantizedBlock = dequantize(quantizedBlock, quantizationTable)
                val idctBlock = upsample(idct(dequantizedBlock), 2, 2)
                for (x in 0 until blockSize * 2) {
                    for (y in 0 until blockSize * 2) {
                        yCbCrDecoded[i * blockSize + x][j * blockSize + y][2] = idctBlock[x][y]
                    }
                }
            }
        }


        // Upsample Cb and Cr channels
        val upsampledCb = upsample(yCbCrDecoded.map { it.map { it[1] }.toIntArray() }.toTypedArray(), 2, 2)
        val upsampledCr = upsample(yCbCrDecoded.map { it.map { it[2] }.toIntArray() }.toTypedArray(), 2, 2)
        for (i in yCbCrDecoded.indices) {
            for (j in yCbCrDecoded[0].indices) {
                yCbCrDecoded[i][j][1] = upsampledCb[i][j]
                yCbCrDecoded[i][j][2] = upsampledCr[i][j]
            }
        }

        // Преобразуем YCbCr обратно в RGB
        val rgbDecoded = yCbCrToRgb(yCbCrDecoded)

        // Создаем изображение и сохраняем его
        val decodedImage = BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
        for (x in 0 until width) {
            for (y in 0 until height) {
                val r = rgbDecoded[x][y][0]
                val g = rgbDecoded[x][y][1]
                val b = rgbDecoded[x][y][2]
                val color = (r shl 16) or (g shl 8) or b
                decodedImage.setRGB(x, y, color)
            }
        }
        ImageIO.write(decodedImage, "jpg", File("D:\\ЛЭТИ\\АиСД\\Второй семестр\\Изображения\\$p.output_$imageName.jpg"))
        p += 10
    //}
}
