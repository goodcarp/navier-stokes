// OCR a PNG with Apple Vision (VNRecognizeTextRequest, accurate mode).
// Usage: swift ocr_vision.swift <file.png> ...
// Output: recognized lines to stdout, one per line, in reading order (top-to-bottom).
import Foundation
import Vision
import AppKit

for path in CommandLine.arguments.dropFirst() {
    guard let img = NSImage(contentsOfFile: path),
          let cg = img.cgImage(forProposedRect: nil, context: nil, hints: nil) else {
        FileHandle.standardError.write("cannot load \(path)\n".data(using:.utf8)!); continue
    }
    let req = VNRecognizeTextRequest()
    req.recognitionLevel = .accurate
    req.usesLanguageCorrection = false
    let handler = VNImageRequestHandler(cgImage: cg, options: [:])
    do { try handler.perform([req]) } catch { FileHandle.standardError.write("fail \(path)\n".data(using:.utf8)!); continue }
    let obs = (req.results ?? []) as [VNRecognizedTextObservation]
    let sorted = obs.sorted { $0.boundingBox.maxY > $1.boundingBox.maxY }
    print("=== \(path) ===")
    for o in sorted { if let c = o.topCandidates(1).first { print(c.string) } }
}
