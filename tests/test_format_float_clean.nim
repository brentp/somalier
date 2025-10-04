import unittest

# Import the function from pedrel module
import ../src/somalierpkg/pedrel

suite "formatFloatClean function":
  test "handles decimal numbers correctly":
    check formatFloatClean(1.234'f32) == "1.234"
    check formatFloatClean(1.2'f32) == "1.2"
    check formatFloatClean(1.1'f32) == "1.1"
    check formatFloatClean(1.0'f32) == "1.0"
    check formatFloatClean(0.5'f32) == "0.5"

  test "handles whole numbers correctly":
    check formatFloatClean(5.0'f32) == "5.0"
    check formatFloatClean(0.0'f32) == "0.0"
    check formatFloatClean(100.0'f32) == "100.0"

  test "handles small decimal values":
    check formatFloatClean(0.1'f32) == "0.1"
    check formatFloatClean(0.01'f32) == "0.01"
    check formatFloatClean(0.001'f32) == "0.001"

  test "handles negative numbers":
    check formatFloatClean(-1.234'f32) == "-1.234"
    check formatFloatClean(-1.2'f32) == "-1.2"
    check formatFloatClean(-1.0'f32) == "-1.0"

  test "handles edge cases with precision":
    check formatFloatClean(1.000001'f32) == "1.000001"
    check formatFloatClean(1.0001'f32) == "1.0001"
    check formatFloatClean(1.000000'f32) == "1.0"

  test "handles very small differences":
    check formatFloatClean(1.0000001'f32) == "1.0"
